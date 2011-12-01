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
use TestUtils;
use Cwd qw/getcwd/;
use File::Basename qw/basename dirname/;
use File::Path qw/make_path remove_tree/;
use File::Spec::Functions qw/canonpath catdir catfile updir/;
use File::Copy qw/copy/;

setup_intermediate_dir();
my $int = $TestUtils::intermediate_dir;

my $testgff = 't/data/test1.gff';
my $copy = catfile($int, basename($testgff));
say "==== $copy";

copy $testgff, $copy;

my ($gene, $chr1, $chr2) = map { catfile($int, "test1-$_.gff") } qw/gene Chr1 Chr2/;

my @a = GFF::Split::split_feature($copy,'gene');
my @b = GFF::Split::split_sequence($copy);

is_deeply(\@a, [ $gene ], "feature");
is_deeply(\@b, [ $chr1, $chr2 ], "sequence");
is_deeply(
    [GFF::Split::split_names($copy, qw/Chr1 Chr2/)], 
    [ $chr1, $chr2 ], 
    "split_names",
);
