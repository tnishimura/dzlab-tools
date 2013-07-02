#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "n|no-index" => \(my $no_index),
);

use File::Basename qw/basename/;
$0 = basename($0);
if (@ARGV == 0){ say STDERR "$0 input.sam"; exit 1}

for my $sam (@ARGV) {
    # my $bam_sorted_prefix = $sam =~ s/\.sam$//r;
    my $bam_sorted_prefix = $sam;
    $bam_sorted_prefix =~ s/\.sam$//;
    my $bam_unsorted      = "$bam_sorted_prefix.unsorted.bam";
    my $bam_sorted        = "$bam_sorted_prefix.bam";
    my $bam_index         = "$bam_sorted_prefix.bam.bai";

    system(qq{ samtools view -bS $sam > $bam_unsorted });
    system(qq{ samtools sort $bam_unsorted $bam_sorted_prefix });
    system(qq{ samtools index $bam_sorted });
    # qq{ samtools tview $bam_sorted }

    unlink $bam_unsorted;
}
