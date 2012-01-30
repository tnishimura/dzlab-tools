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

BEGIN{
    use_ok 'Run::Bowtie';
}

my $reads = "t/data/bs-sequel-test.fastq";
my $dir = setup_intermediate_dir();
my $ref = setup_reference() . '.c2t';
my $output = catfile($dir, "run-bowtie.out");

{
    my ($processed, $aligned, @loglines) = bowtie(
        '-1' => $reads, 
        output => $output, 
        index => $ref,
        verbose => 1,
        seed => 12345,
    );
    is($processed, 5100, 'number processed');
    is($aligned, 3589, 'number aligned');

    my @expected_lines = (
        '# reads processed: 5100',
        '# reads with at least one reported alignment: 3589 (70.37%)',
        '# reads that failed to align: 1511 (29.63%)',
        'Reported 3589 alignments to 1 output stream(s)',
    );

    is_deeply(\@loglines, \@expected_lines, "log lines");
}


