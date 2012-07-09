#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More qw(no_plan);
use Test::Exception;
use TestUtils;
use Cwd qw/getcwd/;
use File::Basename qw/basename dirname/;
use File::Path qw/make_path remove_tree/;
use File::Spec::Functions qw/rel2abs canonpath catdir catfile updir/;
use File::Copy;
use File::Temp qw/tempdir/;
my $tempdir = tempdir(CLEANUP => 1);

use Run::Bowtie;
use Run::BowtieBuild;

my $reads = "t/data/bs-sequel-test.fastq";
my $dir = setup_intermediate_dir();
my $original_ref = setup_reference();
my ($ref) = bowtie_build(file => $original_ref, bs => 'c2t', force => 1);
my $output = catfile($dir, "run-bowtie.out");

{
    my ($processed, $aligned, $suppressed, $reported, @loglines) = bowtie(
        '-1' => $reads, 
        output => $output, 
        index => $ref,
        verbose => 1,
        seed => 12345,
    );
    is($processed, 5100, 'number processed');
    is($aligned, 3589, 'number aligned');
    is($suppressed, 0, 'number suppressed');
    is($reported, 3589, 'total reported');

    my @expected_lines = (
        '# reads processed: 5100',
        '# reads with at least one reported alignment: 3589 (70.37%)',
        '# reads that failed to align: 1511 (29.63%)',
        'Reported 3589 alignments to 1 output stream(s)',
    );

    is_deeply(\@loglines, \@expected_lines, "log lines");
}


{
    my ($processed, $aligned, $suppressed, $reported, @loglines) = bowtie(
        '-1'       => $reads,
        output     => $output,
        index      => $ref,
        verbose    => 1,
        seed       => 12345,
        maxhits    => 10,
        splice     => [5,25],
        readlength => 100,
    );
    is($processed, 5100, 'number processed');
    is($aligned, 4907, 'number aligned');
    is($suppressed, 60, 'number suppressed');
    is($reported, 6137, 'total reported');

    my @expected_lines = (
        '# reads processed: 5100',
        '# reads with at least one reported alignment: 4907 (96.22%)',
        '# reads that failed to align: 133 (2.61%)',
        '# reads with alignments suppressed due to -m: 60 (1.18%)',
        'Reported 6137 alignments to 1 output stream(s)',
    );

    is_deeply(\@loglines, \@expected_lines, "log lines");
}

ok(sub{
    bowtie(
        '-1'       => $reads,
        output     => $output,
        index      => $ref,
        seed       => 12345,
        maxhits    => 10,
        splice     => [5,25],
    )}, "splice without readlength");

dies_ok(sub{
    bowtie(
        '-1'       => $reads,
        output     => $output,
        index      => $ref,
        seed       => 12345,
        maxhits    => 10,
        splice     => [5,25],
        readlength => 100,
        trim5      => 123,
        trim3      => 543,
    )}, "death: splice with trims");

dies_ok(sub{
    bowtie(
        '-1'    => $reads,
        output  => $output,
        index   => $ref,
        seed    => 12345,
        maxhits => 10,
        best    => 1
    )}, "death: max hits with best");
