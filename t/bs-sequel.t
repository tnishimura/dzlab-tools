package bssequel::Test;
use base qw(Test::Class);
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More 'no_plan';
use File::Temp qw/tempdir/;
use File::Basename;
use File::Spec::Functions;
use Launch;
use Carp;

sub setup_files : Test(setup) {
    my $self = shift;
    ++$self->{counter};
    my $gz = 't/TAIR_mini.fas.gz';

    my $test_dir = 't_intermediate';
    my $ref = catfile($test_dir, basename($gz, '.gz'));
    my $counter = $self->{counter};
    if (-e $test_dir && ! -d $test_dir){
        croak "$test_dir is not a directory? dying";
    }
    mkdir $test_dir if ! -d $test_dir;
    launch("gunzip -c $gz > ??", expected => $ref);
    is(8527645, [stat($ref)]->[7], "unzipped reference the correct size");

    my $numreads = 10000;
    my $reads = catfile($test_dir, 'bs-sequel-test.fastq');
    launch("perl ./genome_shear.pl -n $numreads -l 100 -o $reads $ref");

    ok(-f $reads, "reads exist");
   
    $self->{numreads} = $numreads;
    $self->{reads} = $reads;
    $self->{outdir} = tempdir(catfile($test_dir,"bs-sequel-$counter-XXXXX"),CLEANUP => 0);
    $self->{testdir} = $test_dir;
    $self->{ref} = $ref;
}

sub single_ends_runs : Tests(5){
    my $self = shift;
    my ($numreads, $reads, $outdir, $testdir, $reference) = @{$self}{qw/numreads reads outdir testdir ref/};

    launch("./bs-sequel.pl -l $reads -f $reference -b test -d $outdir --no-windowing 2>&1 > /dev/null", verbose => 0);
    ok(-d catfile($outdir, 'single-c'), 'single-c dir exists');
    ok(-d catfile($outdir, 'windows'), 'windows dir exists');
    ok(exists_and_nonempty(catfile($outdir, 'test.gff')), 'test.gff exists');

    my $bowtie = catfile($outdir, "bs-sequel-test.fastq_1-100.eland3.post");
    ok(exists_and_nonempty($bowtie), 'bowtie output exists');
    like(`wc -l $bowtie`, qr/\b$numreads\b/, 'bowtie output has correct number of reads');
}

sub paired_ends_runs : Tests(7){
    my $self = shift;
    my ($numreads, $reads, $outdir, $testdir, $reference) = @{$self}{qw/numreads reads outdir testdir ref/};


    launch("./bs-sequel.pl -ls 1 50 -rs 51 100 -l $reads -f $reference -b test -d $outdir --no-windowing 2>&1 > /dev/null", verbose => 0);
    ok(-d catfile($outdir, 'single-c'), 'single-c dir exists');
    ok(-d catfile($outdir, 'windows'), 'windows dir exists');
    ok(exists_and_nonempty(catfile($outdir, 'test.gff')), 'test.gff exists');

    my $bowtie_left = catfile($outdir, "bs-sequel-test.fastq_1-50.eland3.post");
    my $bowtie_right = catfile($outdir, "bs-sequel-test.fastq_51-100.eland3.post");
    ok(exists_and_nonempty($bowtie_left), 'bowtie_left output exists');
    ok(exists_and_nonempty($bowtie_right), 'bowtie_right output exists');
    like(`wc -l $bowtie_left`, qr/\b$numreads\b/, 'bowtie_left output has correct number of reads');
    like(`wc -l $bowtie_right`, qr/\b$numreads\b/, 'bowtie_right output has correct number of reads');
}

bssequel::Test->runtests;

sub exists_and_nonempty{
    my $file = shift;
    return defined $file &&  -f $file && [stat($file)]->[7];
}
