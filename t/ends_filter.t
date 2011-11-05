#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;

use Test::More qw(no_plan);

my $input = "t/ends_filter.ends";
my $output3 = "t/ends_filter.ends.out.3";
my $output5 = "t/ends_filter.ends.out.5";
my $exons = "t/ends_filter.exon-anno.gff";
my $genes = "t/ends_filter.gene-anno.gff";

system("ends_filter.pl --everywhere -t 100 -5 -b 100 -d 5000 -g $genes -gl ID -e $exons -el Parent -o $output5 $input");
system("ends_filter.pl --everywhere -t 100 -3 -b 100 -d 5000 -g $genes -gl ID -e $exons -el Parent -o $output3 $input");

for my $file ($output3, $output5) {
    open my $fh, '<', $file;
    open my $out, '>', $file . '.filtered';

    my $line1 = <$fh>; chomp $line1;
    my $line2 = <$fh>; chomp $line2;
    my ($at1_label, @at1) = split /\t/, $line1;
    my ($at2_label, @at2) = split /\t/, $line2;

    is($at1_label, "AT1G01010", "line label1");
    is($at2_label, "AT1G01020", "line label2");
    is(scalar(@at1), 100, "number of bins 1");
    is(scalar(@at1), 100, "number of bins 2");
    is_deeply(\@at1, \@at2, "two genes match");

    say $file;
    say $at1_label;
    for (-5 .. 4){ say join "\t", "$_: ", map { shift @at1 } (1 .. 10); }

    say $at2_label;
    for (-5 .. 4){ say join "\t", "$_: ", map { shift @at2 } (1 .. 10); }

    close $fh;
    close $out;
}
