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
    my ($staging, $meta) = prepare_fasta(
        file       => $fasta_file->{file},
        stagingdir => $stagingdir,
    );

    $pm->finish(0, ['fasta', {
                staging => $staging,
                meta    => $meta,
                source  => $fasta_file->{source},
            }]);
}

# GFF lower-density annotation files
for my $gff_file (@{$config->{gff}}) {
    msg($gff_file);
    $pm->start and next;
    my ($staging, @features) = prepare_gff(
        file       => $gff_file->{file},
        stagingdir => $stagingdir,
        source     => $gff_file->{source},
    );
    $pm->finish(0, ['gff', {
            staging => $staging,
            feature => [@features],
            source  => $gff_file->{source},
        }]);
}

# GFF high-density single-c/windows files
for my $gffwig_file (@{$config->{gffwig}}) {
    msg($gffwig_file);
    $pm->start and next;

    my @meta = prepare_gff_to_wig(
        file       => $gffwig_file->{file},
        ctscore    => $gffwig_file->{ctscore},
        source     => $gffwig_file->{source},
        stagingdir => $stagingdir,
        wigdir     => $wigdir,
    );
    $pm->finish(0, ['gffwig', \@meta]);
}

$pm->wait_all_children;

#######################################################################

my %info;
sub collect_staging_files{ # call before calling start()
    my $collected = $_[5];

    die "bug, nothing collected?" unless (ref $collected eq 'ARRAY');
    my ($type, $inforef) = @$collected;

    if (ref $inforef eq 'ARRAY'){
        push @{$info{$type}}, @$inforef;
    }
    else{
        push @{$info{$type}}, $inforef;
    }
};
say Dump(\%info);
DumpFile($config->{output}, \%info);

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
