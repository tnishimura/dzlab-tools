#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use DBI;

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) if ! $opt_database;

my $sth;
my $dbh = DBI->connect("dbi:SQLite:dbname=$opt_database","","", {RaiseError => 1, AutoCommit => 1});

#$dbh->do("create table scores (chr, coord)");
#$dbh->do("create table columns (filename, nickname, num)");

print "chromosome\tstart\tend";
for my $colspec (@{$dbh->selectall_arrayref("select nickname from columns order by num")}){
    my ($nickname) = @$colspec;
    print "\t${nickname}_score";
    print "\t${nickname}_ct";
}
print "\n";


my $select = $dbh->prepare("select * from scores order by chr, coord");
$select->execute();
while (defined(my $row_ref = $select->fetchrow_arrayref())){
    my @row = @$row_ref;
    if (@row < 4 || @row % 2 != 0){
        die "something wrong, unexpected number of col in $opt_database.scores";
    }
    my $chr = shift @row;
    my $coord = shift @row;
    print "$chr\t$coord\t$coord";
    while (scalar @row){
        my $c = shift @row;
        my $t = shift @row;
        $c //= 0;
        $t //= 0;
        my $total = $c + $t;
        if ($total == 0){
            print "\t.\t.";
        }
        else{
            my $score = $c/$total;
            printf "\t%.7f\t%d", $score, $total;
        }
    }
    print "\n";
}
$dbh->disconnect;

=head1 NAME
 
single_c_mashup_create.pl
 
=head1 SYNOPSIS

 single_c_mashup_create.pl

=head1 OPTIONS

=over

=item  -d <file> | --database <file>

=item  -o <file> | --output <file>

=for Euclid
    file.default:     '-'

=back

=cut

