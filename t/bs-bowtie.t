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
my $num_junk = 1000;
my $total = $num_reads + $num_junk;
my $bowtieout = catfile($test_dir, "bs-bowtie-test");
my $testfastq = catfile($test_dir, "bs-bowtie-sheared.fastq");

launch("perl ./genome_shear.pl -j $num_junk -n $num_reads -l 100 -o $testfastq $ref");
launch("perl ./bs-bowtie -f $ref -r $testfastq -o $bowtieout");
like(`wc -l $bowtieout`, qr/$total/, "correct number of reads output from bs-bowtie");

my $NMcount = 0;
my %mismatches;
{
    open my $fh, '<', $bowtieout;
    while (defined(my $line = <$fh>)){
        chomp $line;
        my @split = split /\t/, $line;
        if ($split[2] eq 'NM'){
            $NMcount++;
        }
    }
    close $fh;
}

ok($NMcount >= $num_junk, "bowtie output has at least as many junks as inputed");

