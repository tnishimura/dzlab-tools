#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;

END {close STDOUT}
$| = 1;

while (defined(my $line = <ARGV>)){
    chomp $line;
    my @split = split /\t/, $line;

    if ($split[8] =~ s/;n=(\d+)//){
        ($split[5], $split[8]) = ($1, $split[8] . ";n=" . int($split[5]));
    }
    say join "\t", @split;
}
__DATA__
chr1    .       w1      -66     -66     0.0000  .       .       c=0;t=0;n=3
chr1    .       w1      -65     -65     0.0000  .       .       c=0;t=0;n=2
chr1    .       w1      -60     -60     0.0000  .       .       c=0;t=0;n=2
chr1    .       w1      -57     -57     0.0000  .       .       c=0;t=0;n=1
chr1    .       w1      -52     -52     0.0000  .       .       c=0;t=0;n=1
chr1    .       w1      -51     -51     0.0000  .       .       c=0;t=0;n=5
chr1    .       w1      -48     -48     0.0000  .       .       c=0;t=0;n=2
chr1    .       w1      -47     -47     0.0000  .       .       c=0;t=0;n=1
chr1    .       w1      -44     -44     0.0000  .       .       c=0;t=0;n=3
chr1    .       w1      -35     -35     0.0000  .       .       c=0;t=0;n=7
chr1    .       w1      -34     -34     0.0000  .       .       c=0;t=0;n=4
chr1    .       w1      -26     -26     0.0000  .       .       c=0;t=0;n=3
chr1    .       w1      -19     -19     0.0000  .       .       c=0;t=0;n=5
chr1    .       w1      -11     -11     0.0000  .       .       c=0;t=0;n=1
chr1    .       w1      -7      -7      0.0000  .       .       c=0;t=0;n=1
chr1    .       w1      -6      -6      0.0000  .       .       c=0;t=0;n=1
chr1    .       w1      -5      -5      0.0000  .       .       c=0;t=0;n=2
chr1    .       w1      -4      -4      0.0000  .       .       c=0;t=0;n=8

