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

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(exists_and_nonempty setup_reference);

sub setup_reference{
    my $test_dir = shift;
    my $gz = 't/TAIR_mini.fas.gz';
    my $ref = catfile($test_dir, basename($gz, '.gz'));
    if (-e $test_dir && ! -d $test_dir){
        croak "$test_dir is not a directory? dying";
    }
    mkdir $test_dir if ! -d $test_dir;
    launch("gunzip -c $gz > ??", expected => $ref, verbose => 0);
    launch("./bs-bowtie-build -c2t $ref", verbose => 0);
    launch("./bs-bowtie-build -g2a $ref", verbose => 0);
    is(8527645, [stat($ref)]->[7], "unzipped reference the correct size");
    return $ref;
}

sub exists_and_nonempty{
    my $file = shift;
    return defined $file &&  -f $file && [stat($file)]->[7];
}

1;

