#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/../lib";
use GFF::Split;
use Test::More qw(no_plan);

my @a = GFF::Split::split_feature('t/test1.gff');
my @b = GFF::Split::split_sequence('t/test1.gff');

is_deeply(\@a, [ 't/test1-gene.gff' ], "feature");
is_deeply(\@b, [ 't/test1-CHR1.gff', 't/test1-CHR2.gff' ], "sequence");
is_deeply(
    [GFF::Split::split_names('t/test1.gff', qw/CHR1 CHR2/)], 
    [ 't/test1-CHR1.gff', 't/test1-CHR2.gff' ], 
    "split_names",
);
