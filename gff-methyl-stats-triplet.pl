#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use Getopt::Long;

use lib "$FindBin::Bin/lib";
use GFF::Statistics qw/methylation_stats/;
use YAML qw/Dump/;

my $result = GetOptions (
    "wp=s" => \(my $wanted_percentiles = ".05,.25,.50,.75,.95"),
);
my @wp = split /,/, $wanted_percentiles;
usage() if (@ARGV != 3 || ! $result);
my %files = (
    cg => $ARGV[0],
    chg => $ARGV[1],
    chh => $ARGV[2],
);

my %total_nuc_ct;
my %total_chr_ct;
my %total_mit_ct;

use Tie::IxHash;
tie my %total_stats, "Tie::IxHash", qw/cg chg chh/;

while (my ($name,$file) = each %files) {
    my ($stats, $nuc_ct, $chr_ct, $mit_ct)  = methylation_stats($file);
    $stats->{file} = $file;

    $total_stats{$name} = $stats;
    while (my ($bin,$count) = each %$nuc_ct) {
        $total_nuc_ct{$bin} += $count;
    }
    while (my ($bin,$count) = each %$chr_ct) {
        $total_chr_ct{$bin} += $count;
    }
    while (my ($bin,$count) = each %$mit_ct) {
        $total_mit_ct{$bin} += $count;
    }
}

$total_stats{'total'} = {
    nuc_ct_percentile => GFF::Statistics::histpercentiles(\%total_nuc_ct, @wp),
    chr_ct_percentile => GFF::Statistics::histpercentiles(\%total_chr_ct, @wp),
    mit_ct_percentile => GFF::Statistics::histpercentiles(\%total_mit_ct, @wp),
    nuc_ct_mean   => GFF::Statistics::histmean(\%total_nuc_ct),
    chr_ct_mean   => GFF::Statistics::histmean(\%total_chr_ct),
    mit_ct_mean   => GFF::Statistics::histmean(\%total_mit_ct),
};

say Dump \%total_stats;

sub usage{
    say STDERR "$0 cg chg chh";
    exit 1;
}

