package TestUtils;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Launch;
use File::Basename;
use File::Spec::Functions;
use Test::More;
use File::Path qw/make_path/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(setup_intermediate_dir exists_and_nonempty setup_reference test_gff);
our $intermediate_dir = 't_intermediate';

sub setup_intermediate_dir{
    make_path $intermediate_dir;
    return $intermediate_dir;
}

sub setup_reference{
    my $test_dir = shift // $intermediate_dir;
    my $do_methylation = shift;

    my $gz = 't/data/TAIR_mini.fas.gz';
    my $ref = catfile($test_dir, basename($gz, '.gz'));
    if (-e $test_dir && ! -d $test_dir){
        croak "$test_dir is not a directory? dying";
    }
    mkdir $test_dir if ! -d $test_dir;
    launch("gunzip -c $gz > ??", expected => $ref);

    if ($do_methylation){
        launch("./bs-bowtie-build -c2t $ref");
        launch("./bs-bowtie-build -g2a $ref");
    }

    is(8527645, [stat($ref)]->[7], "unzipped reference the correct size");
    return $ref;
}

sub exists_and_nonempty{
    my $file = shift;
    return defined $file &&  -f $file && [stat($file)]->[7];
}

sub test_gff{
    my $line_count = shift // 1_000_000;
    my $file = shift // catfile(setup_intermediate_dir(), 'test.gff');
    if (-s $file){
        return $file;
    }

    open my $fh, '>', $file ;
    for (1 .. $line_count) {
        my $chr   = 'chr' . int(rand(100));
        my $start = int(rand(10_000_000));
        my $end   = $start + int(rand(10_000));
        my $c     = int(rand(1000));
        my $t     = int(rand(1000));
        my $n     = int(rand(1000));
        say $fh join "\t", 
            $chr,
            q{.},
            q{.},
            $start, 
            $end,
            rand(1000),
            (rand() > .5 ? q{+} : q{-}),
            q{.},
            "c=$c;t=$t;n=$n"
        ;
    }
    close $fh;

    return $file;
}

1;

