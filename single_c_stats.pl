#!/usr/bin/env perl
use 5.10.0;
use strict;
use warnings FATAL => "all";
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use List::Util qw/sum/;

END {close STDOUT}
$|=1;


say "file\tcount\tmean\tmedian";


for my $file (@ARGV) {
    my $p = GFF::Parser->new(file => $file);
    my @cplust;
    while (defined(my $gff = $p->next())){
        my ($c, $t) = ($gff->get_column('c'), $gff->get_column('t'));
        next unless defined $c && defined $t;
        push @cplust, $c+$t;
    }
    @cplust = sort { $a <=> $b } @cplust;
    my $length = scalar @cplust;
    my $mean   = $length ? sum(@cplust) / $length : 0;
    my $median = $length ? $cplust[int($length/2)]: 0;
    printf "%s\t%d\t%.4f\t%.4f\n", $file, $length, $mean, $median;
    @cplust = ();
}




