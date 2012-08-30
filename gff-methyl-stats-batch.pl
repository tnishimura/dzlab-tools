#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use IO::All;
use Parallel::ForkManager;

use Pod::Usage;
use Getopt::Long;

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "p|parallel=i" => \(my $p),
);
if (! $result || (! @ARGV && -t STDIN)){
    pod2usage(-verbose => 2, -noperldoc => 1);
}

my $pm = Parallel::ForkManager->new($p);

for my $dir (@ARGV) {
    for my $infile (glob("$dir/single-c/*all*.gff")) {
        $pm->start and next;
        pdo{
            my $methfile = $infile =~ s/gff$/methstats.txt/r;
            my $cmd = "gff-methyl-stats.pl -d mstats $infile > $methfile";
            say $cmd;
            system($cmd);
        }
        $pm->finish; 
    }
}
$pm->wait_all_children;

=head1 gff-methyl-stats-batch.pl 

Run after ends_analysis_batch.pl

Usage examples:

 gff-methyl-stats-batch.pl -p 8

=cut

