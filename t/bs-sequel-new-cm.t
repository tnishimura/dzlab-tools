package bssequel::Test;
use base qw(Test::Class);
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More tests => 36;
#use Test::More 'no_plan';
use File::Temp qw/tempdir/;
use File::Basename;
use File::Spec::Functions;
use Launch;
use Carp;
use TestUtils;
use YAML qw/LoadFile/;
use List::MoreUtils qw/all/;

$Launch::VERBOSE = $ENV{HARNESS_VERBOSE};

sub setup_files : Test(setup => 3) {
    my $self = shift;

    my $test_dir = 't_intermediate';
    my $ref = setup_reference($test_dir);

    # just testing that it runs, so few reads
    my $numreads = 50000;
    my $junk = 1000;
    my $reads = catfile($test_dir, 'bs-sequel-test.fastq');

    ok( launch("./genome_shear.pl -j $junk -n $numreads -l 100 -o $reads $ref"),
            "genome_shear.pl ran" );
    ok(-f $reads, "reads exist");
   
    $self->{numreads} = $numreads + $junk;
    $self->{reads} = $reads;
    $self->{methsites} = catfile($test_dir, 'bs-sequel-test.fastq.methsites.gff');
    $self->{testdir} = $test_dir;
    $self->{ref} = $ref;
}

sub single_ends_runs : Tests(14){
    my $self = shift;
    my ($numreads, $reads, $testdir, $reference) = @{$self}{qw/numreads reads testdir ref/};
    my $outdir = tempdir(catfile($testdir,"bs-sequel-single-XXXXX"),CLEANUP => 0);
    my $bowtie = catfile($outdir, "bs-sequel-test.fastq_1-100.eland3.post");
    my $correlate = catfile($outdir, "test.gff");

    # run bs-sequel
    ok(launch("./bs-sequel.pl --new-cm -l $reads -f $reference -b test -d $outdir --no-windowing 2>&1 > /dev/null"),
        "bs-sequel.pl survived");

    # basics
    ok(-d catfile($outdir, 'single-c'), 'single-c dir exists');
    ok(-d catfile($outdir, 'windows'), 'windows dir exists');
    ok(exists_and_nonempty(catfile($outdir, 'test.gff')), 'test.gff exists');

    # bowtie file
    ok(exists_and_nonempty($bowtie), 'bowtie output exists');
    like(`wc -l $bowtie`, qr/\b$numreads\b/, 'bowtie output has correct number of reads');

    # single c
    my @singlec = glob(catfile($outdir, 'single-c', '*merged'));
    is(scalar @singlec, 3, '3 single-c files');
    for my $sc (@singlec) {
        $self->single_c_confirm($sc);
    }
    $self->single_c_diff(@singlec);

    # correlated single ends output
    ok(launch("./correlate_check.pl -r $reference $correlate"), "correlate_check.pl passed");
}

sub paired_ends_runs : Tests(16){
    my $self = shift;
    my ($numreads, $reads, $testdir, $reference) = @{$self}{qw/numreads reads testdir ref/};
    my $outdir = tempdir(catfile($testdir,"bs-sequel-paired-XXXXX"),CLEANUP => 0);
    my $bowtie_left = catfile($outdir, "bs-sequel-test.fastq_1-60.eland3.post");
    my $bowtie_right = catfile($outdir, "bs-sequel-test.fastq_61-100.eland3.post");
    my $correlate = catfile($outdir, "test.gff");

    # run bs-sequel
    ok(launch("./bs-sequel.pl --new-cm -ls 1 60 -rs 61 100 -l $reads -f $reference -b test -d $outdir --no-windowing 2>&1 > /dev/null"),
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

    # single c
    my @singlec = glob(catfile($outdir, 'single-c', '*merged'));
    is(scalar @singlec, 3, '3 single-c files');
    for my $sc (@singlec) {
        $self->single_c_confirm($sc);
    }
    $self->single_c_diff(@singlec);

    # correlated single ends output
    ok(launch("./correlate_check.pl -r $reference $correlate"), "correlate_check.pl passed");
}

sub single_c_confirm{ # 3
    my ($self, $singlec) = @_;
    ok(launch("./single_c_confirm.pl -r $self->{ref} $singlec"),
        "single-c contains only C positions");
}

sub single_c_diff{ # 3
    my ($self, @singlecs) = @_;
    my $diff = catfile($self->{testdir}, 'scdiff');
    ok(launch("./single_c_diff.pl $self->{methsites} @singlecs > $diff"), 'single_c_diff.pl ran');
    my $output = LoadFile($diff);
    say Dumper $output;
    my $has_keys = all {exists $output->{$_}} qw/ correct error found_and_underscore non_methylated_detected overscore total_detected_c total_real_c/;
    ok($has_keys, 'has all keys');
    my ($detected, $real) = @{$output}{qw/ total_detected_c total_real_c/};
    ok($detected/$real > .9, 'at least 90% detected');
}

bssequel::Test->runtests;

