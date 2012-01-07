#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Perl6::Slurp;
use Test::More qw(no_plan);
$|=1;

BEGIN {use_ok "DZUtil", 'downsample';}

my $data = 't/data/downsample_input.txt';
my $multiple = 3;

my ($tempfile, $sampled, $total) = downsample $data, .5, $multiple;
ok($sampled < $total, "sampled < total");

my @lines = split /\n/, slurp($tempfile);
ok(@lines == $sampled, "correct number of sampled lines"); 
is(0, @lines % $multiple, "contents multiple of 3");

for (my $i = 0 ; $i <= $#lines ; $i+=$multiple){
    my ($x,$y,$z) = @lines[$i .. $i + $multiple - 1];
    ok($x eq $y && $y eq $z, "triplets are equal");
}
