package bssequel::Test;
use base qw(Test::Class);
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More tests => 22;
use File::Temp qw/tempdir/;
use File::Basename;
use File::Spec::Functions;
use Launch;
use Carp;
use TestUtils;

sub setup_files : Test(setup => 3) {
    my $self = shift;

    ++$self->{counter};
    my $counter = $self->{counter};

    my $test_dir = 't_intermediate';
    my $ref = setup_reference($test_dir);

    # just testing that it runs, so few reads
    my $numreads = 1000;
    my $junk = 100;
    my $reads = catfile($test_dir, 'bs-sequel-test.fastq');

    ok( launch("./genome_shear.pl -j $junk -n $numreads -l 100 -o $reads $ref"),
            "genome_shear.pl ran" );
    ok(-f $reads, "reads exist");
   
    $self->{numreads} = $numreads + $junk;
    $self->{reads} = $reads;
    $self->{outdir} = tempdir(catfile($test_dir,"bs-sequel-$counter-XXXXX"),CLEANUP => 0);
    $self->{testdir} = $test_dir;
    $self->{ref} = $ref;
}

sub single_ends_runs : Tests(7){
    my $self = shift;
    my ($numreads, $reads, $outdir, $testdir, $reference) = @{$self}{qw/numreads reads outdir testdir ref/};
    my $bowtie = catfile($outdir, "bs-sequel-test.fastq_1-100.eland3.post");
    my $correlate = catfile($outdir, "test.gff");

    # run bs-sequel
    ok(launch("./bs-sequel.pl -l $reads -f $reference -b test -d $outdir --no-windowing 2>&1 > /dev/null", verbose => 0),
        "bs-sequel.pl survived");

    # basics
    ok(-d catfile($outdir, 'single-c'), 'single-c dir exists');
    ok(-d catfile($outdir, 'windows'), 'windows dir exists');
    ok(exists_and_nonempty(catfile($outdir, 'test.gff')), 'test.gff exists');

    # bowtie file
    ok(exists_and_nonempty($bowtie), 'bowtie output exists');
    like(`wc -l $bowtie`, qr/\b$numreads\b/, 'bowtie output has correct number of reads');

    # correlated single ends output
    ok(launch("./correlate_check.pl -r $reference $correlate"), "correlate_check.pl passed");
}

sub paired_ends_runs : Tests(9){
    my $self = shift;
    my ($numreads, $reads, $outdir, $testdir, $reference) = @{$self}{qw/numreads reads outdir testdir ref/};
    my $bowtie_left = catfile($outdir, "bs-sequel-test.fastq_1-60.eland3.post");
    my $bowtie_right = catfile($outdir, "bs-sequel-test.fastq_61-100.eland3.post");
    my $correlate = catfile($outdir, "test.gff");

    # run bs-sequel
    ok(launch("./bs-sequel.pl -ls 1 60 -rs 61 100 -l $reads -f $reference -b test -d $outdir --no-windowing 2>&1 > /dev/null", verbose => 0),
        "bs-sequel.pl survived");

    # basics
    ok(-d catfile($outdir, 'single-c'), 'single-c dir exists');
    ok(-d catfile($outdir, 'windows'), 'windows dir exists');
    ok(exists_and_nonempty(catfile($outdir, 'test.gff')), 'test.gff exists');

    # bowtie
    ok(exists_and_nonempty($bowtie_left), 'bowtie_left output exists');
    ok(exists_and_nonempty($bowtie_right), 'bowtie_right output exists');
    like(`wc -l $bowtie_left`, qr/\b$numreads\b/, 'bowtie_left output has correct number of reads');
    like(`wc -l $bowtie_right`, qr/\b$numreads\b/, 'bowtie_right output has correct number of reads');

    # correlated single ends output
    ok(launch("./correlate_check.pl -r $reference $correlate"), "correlate_check.pl passed");
}

bssequel::Test->runtests;

