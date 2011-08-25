#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Launch;
use File::Basename;
use File::Spec::Functions;
use Test::More qw(no_plan);

# setup. probably belongs in a function
my $test_dir = 't_intermediate';
my $gz = 't/TAIR_mini.fas.gz';
my $ref = catfile($test_dir, basename($gz, '.gz'));
mkdir $test_dir if ! -d $test_dir;
launch("gunzip -c $gz > ??", expected => $ref);


my $num_reads = 10000;
my $bowtieout = catfile($test_dir, "bs-bowtie-test");
my $testfastq = catfile($test_dir, "bs-bowtie-sheared.fastq");

launch("perl ./genome_shear.pl -n $num_reads -l 100 -o $testfastq $ref");
launch("perl ./bs-bowtie -f $ref -r $testfastq -o $bowtieout");
like(`wc -l $bowtieout`, qr/$num_reads/, "correct number of redas output from bs-bowtie");

