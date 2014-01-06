#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use IO::File;
use IO::Handle;
use File::Temp;
use File::Copy qw/move/;
use DBI;
use Pod::Usage;
use Getopt::Long;
use List::Util qw/first max min shuffle sum/;
use File::Temp;

use FindBin;
use lib "$FindBin::Bin/lib";
use FastaReader;

my $result = GetOptions (
    "reference|r=s" => \(my $reference),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || ! $reference);  

my ($dbh, $select) = make_index($reference);

while (defined(my $line = <ARGV>)){
    chomp $line;
    my @f = split /\t/, $line;
    next if @f != 9;
    my ($seq, $start) = @f[0,3];
    $f[6] = get_strand($select, uc($seq), $start);
    say join "\t", @f;
}

sub make_index{
    my $reference = shift;

    my $sqlite_file = "$reference.methstrand.sqlite3";
    if (! -f $sqlite_file){
        my $tmpfh = File::Temp->new( TEMPLATE => "$sqlite_file.tmpXXXXX", UNLINK => 0);
        my $tmp_filename = $tmpfh->filename;

        my $fr = FastaReader->new(file => $reference, slurp => 1,);
        my %seqlen = $fr->sequence_lengths;

        my $dbh = DBI->connect("dbi:SQLite:dbname=$tmp_filename","","", {RaiseError => 1, AutoCommit => 0});
        $dbh->do(q{create table tab (seq, pos integer, strand)});
        $dbh->do(q{create index idx1 on tab (seq, pos)});
        my $insert_sth = $dbh->prepare(q{insert into tab (seq, pos, strand) values (?,?,?)});

        my $step = 1_000_000;
        for my $seqid (sort keys %seqlen) {
            my $max = $seqlen{$seqid};
            $seqid = uc $seqid;
            my $start = 1;
            while (){
                if ($start > $max){
                    last;
                }
                my $end = min($start + $step - 1, $max);;

                my @arr = split //, $fr->get($seqid, $start, $end);
                if (@arr != $end - $start + 1){
                    die "\$fr->get($seqid, $start, $end) not expected size";
                }
                for (my $i = $start ; $i <= $end; ++$i){
                    my $base = uc $arr[$i - $start];
                    if ($base eq 'C'){
                        $insert_sth->execute($seqid, $i, '+');
                    }
                    elsif ($base eq 'G'){
                        $insert_sth->execute($seqid, $i, '-');
                    }
                }

                $start += $step - 1;
            }
            $dbh->commit;
        }
        $dbh->commit;
        $dbh->disconnect;

        move $tmp_filename, $sqlite_file;
    }

    my $dbh = DBI->connect("dbi:SQLite:dbname=$sqlite_file","","", {RaiseError => 1, ReadOnly => 1,});
    my $select = $dbh->prepare(q{select strand from tab where seq = ? and pos = ?});
    return ($dbh, $select);
}

sub get_strand{
    my $sth = shift;
    $sth->execute(@_);
    if (my $row = $sth->fetchrow_hashref) {
        return $row->{strand};
    }
    return '.';
}

=head1 gff-add-strand.pl 

Usage examples:

 gff-add-strand.pl -r TAIR_reference.fas input.gff

=cut


