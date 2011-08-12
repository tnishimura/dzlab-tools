#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;

END {close STDOUT}
$| = 1;

if (2 != @ARGV){
    say "$0 - subtract X's scores from Y's and output new GFF";
    say "usage: $0 X.gff Y.gff > result.gff";
    exit 1;
}

my @files = @ARGV;

open my $left_fh, '<', $files[0];
open my $right_fh, '<', $files[1];

while (
    defined (my $left_line = <$left_fh>) and 
    defined (my $right_line = <$right_fh>) 
){
    chomp($left_line, $right_line);
    my @left = split /\t/, $left_line;
    my @right = split /\t/, $right_line;

    say join "\t", @left[0 .. 4], (($left[5] // 0) - ($right[5] //0)), @right[6..8];
}
