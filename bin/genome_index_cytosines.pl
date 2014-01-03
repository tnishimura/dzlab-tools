#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/../lib";
use FastaReader;
use Pod::Usage;
use Getopt::Long;
use DBI;

END {close STDOUT}
$| = 1;

my $result = GetOptions (
    "reference|r=s" => \(my $reference),
    "output|o=s"    => \(my $output),
    "force|f"       => \(my $force),
);
if (!$result || !$reference){
    say "usage: genome_index_cytosines.pl ...";
    exit 1;
}

$output //= $reference . ".cindex";

{
    if (-f $output){
        if ($force){
            say "$output already existed, overwriting.";
            unlink $output;
        }
        else{
            say "$output already exists, not overwriting unless -f is given. exitting.";
            exit 1;
        }
    }
    my $dbh = DBI->connect("dbi:SQLite:dbname=$output","","", {RaiseError => 1, AutoCommit => 1});
    $dbh->do("PRAGMA automatic_index = OFF");
    $dbh->do("PRAGMA journal_mode = OFF");
    $dbh->do("PRAGMA cache_size = 80000");
    $dbh->do("create table methyl (seq, coord integer, context, c integer default 0, t integer default 0)");
    $dbh->{AutoCommit} = 0;
    my $sth = $dbh->prepare("insert into methyl (seq,coord,context) values (?,?,?)");

    my $counter=0;
    sub add{
        $sth->execute(@_);
        if ($counter++ % 20000 == 0){
            $dbh->commit();
            say $counter;
        }
    }

    sub finalize{
        $dbh->commit();
        $dbh->do("create index idx on methyl (seq,coord)");
    }
}

my $fr = FastaReader->new(file => $reference, slurp => 1);
my $regex = qr/(C|G)/;

for my $seqname ($fr->sequence_list()) {
    my $sequence = $fr->get($seqname, undef, undef, rc => 0, coord => 'f');
    while ($sequence =~ /$regex/g){
        my $position = $-[0] + 1;
        add($seqname,
            $position,
            $fr->get_context($seqname, $position, base => 1 =>  rc => $1 eq 'G'));
    }
}
finalize();
