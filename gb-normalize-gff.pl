#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Scalar::Util qw/looks_like_number/;

END {close STDOUT}
$| = 1;

while (defined(my $line = <ARGV>)){
    $line =~ tr/\n\r//d;
    my @F = split /\t/, $line, 9;
    if (@F ==9 and looks_like_number($F[3]) && looks_like_number($F[4]) && $F[3] <= $F[4]){
        $F[0] = lc $F[0]; say join "\t", @F
    }
}
