#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use List::Util qw/sum/;
use Statistics::Descriptive;

END {close STDOUT}
$| = 1;


for my $file (@ARGV) {
    my $count = 0;
    my %seqs;
    my %features;

    my $lengths = Statistics::Descriptive::Full->new();

    my $parser = GFF::Parser->new(file => $file, normalize => 0);
    while (defined(my $gff = $parser->next())){
        my ($seq, $feature) = ($gff->sequence(), $gff->feature());
        $seqs{$seq}++;
        $features{$feature}++;
        $count++;
        $lengths->add_data($gff->end - $gff->start + 1);
    }
    say "### $file";
    printf("# Length mean std (min, 25%%, 50%%, 75%%, max): %.1f %.1f (%d,%d)\n", 
        $lengths->mean(), 
        $lengths->standard_deviation(),
        $lengths->min(), 
        $lengths->max(), 
    );
    printf("# Sequences (%d total):\n", scalar(keys %seqs));
    for my $seq (sort keys %seqs) { say "$seq\t$seqs{$seq}"; }

    printf("# Features (%d total):\n", scalar(keys %features));
    for my $fet (sort keys %features) { say "$fet\t$features{$fet}"; }
}

