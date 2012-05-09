#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More qw(no_plan);
use Test::Exception;

use MethylCounter;
use FastaReader;
use TestUtils;
use Cwd qw/getcwd/;
use File::Basename qw/basename dirname/;
use File::Path qw/make_path remove_tree/;
use File::Spec::Functions qw/rel2abs canonpath catdir catfile updir/;
use File::Copy;
use File::Temp qw/tempdir/;

my $tmpdir = setup_intermediate_dir();

my $ref = FastaReader->new(file => setup_reference($tmpdir, 0), slurp => 1);

my $mc = MethylCounter->new(reference_genome => $ref, output_prefix => catfile($tmpdir, "meow"));

my $base = 1;

isa_ok($mc, 'MethylCounter');

$mc->record_methylation('chr1', 1000, 'C');
$mc->commit();

my $pdl = $mc->{bigarrays}{'CHR1'}{C}{pdl};

is(1, $pdl->at(1000 - $base), "single insertion point");
is(1, $pdl->sum(), "single insertion sum");

$mc->record_methylation('chr1', 1000, 'C');
$mc->record_methylation('chr1', 2000, 'T');
$mc->record_methylation('chr1', 3000, 'C');
$mc->commit();

is(2, $pdl->at(1000 - $base), "single insertion point again");

