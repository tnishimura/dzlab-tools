#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use YAML qw/LoadFile DumpFile Dump/;
use Parallel::ForkManager;
use List::MoreUtils qw/notall/;
use File::Path qw/make_path/;

use FindBin;
use lib "$FindBin::Bin/lib";
use GBUtil;
use GBUtil::InputFile::MethylGFF;
use GBUtil::InputFile::GFF;
use GBUtil::InputFile::Fasta;

my $config_file = shift // usage();
my $config = LoadFile($config_file);
my $pm = Parallel::ForkManager->new($config->{parallel});

my $stagingdir = $config->{stagingdir};
my $wigdir = $config->{wigdir};

if (notall { defined $_ } qw/stagingdir wigdir/){
    usage();
}

make_path($stagingdir);
make_path($wigdir);

# collect files from 
$pm->run_on_finish(\&collect_staging_files);

# FASTA files
for my $fasta_file (@{$config->{fasta}}) {
    msg($fasta_file);
    $pm->start and next;
    my $fasta_track = GBUtil::InputFile::Fasta->new(
        file        => $fasta_file->{file},
        staging_dir => $stagingdir,
        source      => $fasta_file->{source}
    );
    $fasta_track->convert;
    $pm->finish(0, $fasta_track);
}

# GFF lower-density annotation files
for my $gff_file (@{$config->{gff}}) {
    msg($gff_file);
    $pm->start and next;
    my $gff_track = GBUtil::InputFile::GFF->new(
        file        => $gff_file->{file},
        staging_dir => $stagingdir,
        source      => $gff_file->{source},
    );
    $gff_track->convert();

    $pm->finish(0, $gff_track);
}

# GFF high-density single-c/windows files
for my $methylgff_file (@{$config->{methylgff}}) {
    msg($methylgff_file);
    $pm->start and next;

    my $meta_track = GBUtil::InputFile::MethylGFF->new(
        file        => $methylgff_file->{file},
        staging_dir => $stagingdir,
        source      => $methylgff_file->{source},

        wig_dir     => $wigdir,
        ctscore     => $methylgff_file->{ctscore},
    );
    $meta_track->convert();

    $pm->finish(0, $meta_track);
}

$pm->wait_all_children;

#######################################################################
# collect from forks

my @tracks;
sub collect_staging_files{ # call before calling start()
    my $collected = $_[5];
    die "bug, nothing collected?" unless (defined $collected);
    push @tracks, $collected;
};

#######################################################################
# output

say Dump(\@tracks);
DumpFile($config->{stagingconf}, \@tracks);

#######################################################################

sub usage{
    say STDERR "gb-prepare.pl config.yaml";
    exit 1;
}

sub msg{
    say STDERR "[+] " . Dumper shift;
}

__DATA__
# example

parallel:    4
stagingdir:  /home/toshiro/demeter/staging
wigdir:      /home/toshiro/demeter/wig

output:      /home/toshiro/demeter/meow.yaml

fasta:
  - file:    /home/toshiro/genomes/AT/TAIR_reference.fas
    source:  TAIR8
gff:
  - file:    /home/toshiro/annotations/AT/gmod/TAIR8_gmod.gff
    source:  TAIR8
methylgff:
  - file:    /home/toshiro/GEO-Submission-AT-Demeter-2012/windows/at-endosperm-ler_fie_x_col_wt/windows-Col/all.cg-col.w50.gff
    source:  at-en-lerfie-x-col-wt-cg
  - file:    /home/toshiro/GEO-Submission-AT-Demeter-2012/windows/at-endosperm-ler_fie_x_col_wt/windows-Col/all.chg-col.w50.gff
    source:  at-en-lerfie-x-col-wt-chg
  - file:    /home/toshiro/GEO-Submission-AT-Demeter-2012/windows/at-endosperm-ler_fie_x_col_wt/windows-Col/all.chh-col.w50.gff
    source:  at-en-lerfie-x-col-wt-chh
