#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use YAML qw/LoadFile/;
use Parallel::ForkManager;
use List::MoreUtils qw/notall/;
use File::Path qw/make_path/;
use FindBin;
use lib "$FindBin::Bin/lib";
use GBUtil;

my $config_file = shift // usage();
my $config = LoadFile($config_file);
my $pm = Parallel::ForkManager->new($config->{parallel});

my $stagingdir = $config->{stagingdir};
my $wigdir = $config->{wigdir};

if (notall { defined $_ } qw/stagingdir wigdir/){
    usage();
}

# FASTA files
for my $fasta_file (@{$config->{fasta}}) {
    msg($fasta_file);
    $pm->start and next;
    prepare_fasta(
        file     => $fasta_file->{file},
        stagingdir => $stagingdir,
    );
    $pm->finish; 
}

# GFF lower-density annotation files
for my $gff_file (@{$config->{gff}}) {
    msg($gff_file);
    $pm->start and next;
    prepare_gff(
        file     => $gff_file->{file},
        stagingdir => $stagingdir,
    );
    $pm->finish; 
}

# GFF high-density single-c/windows files
for my $gffwig_file (@{$config->{gffwig}}) {
    msg($gffwig_file);
    $pm->start and next;
    prepare_gff_to_wig(
        file     => $gffwig_file->{file},
        ctscore  => $gffwig_file->{ctscore},
        source   => $gffwig_file->{source},
        stagingdir => $stagingdir,
        wigdir   => $wigdir,
    );
    $pm->finish; 
}

$pm->wait_all_children;

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
fasta:
  - file:    /home/toshiro/genomes/AT/TAIR_reference.fas
gff:
  - file:    /home/toshiro/annotations/AT/gmod/TAIR8_gmod.gff
gffwig:
  - file:    /home/toshiro/GEO-Submission-AT-Demeter-2012/windows/at-endosperm-ler_fie_x_col_wt/windows-Col/all.cg-col.w50.gff
    source:  at-en-lerfie-x-col-wt-cg
    ctscore: 0 
  - file:    /home/toshiro/GEO-Submission-AT-Demeter-2012/windows/at-endosperm-ler_fie_x_col_wt/windows-Col/all.chg-col.w50.gff
    source:  at-en-lerfie-x-col-wt-chg
    ctscore: 0 
  - file:    /home/toshiro/GEO-Submission-AT-Demeter-2012/windows/at-endosperm-ler_fie_x_col_wt/windows-Col/all.chh-col.w50.gff
    source:  at-en-lerfie-x-col-wt-chh
    ctscore: 0 
