#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Pod::Usage;
use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser::Correlated;
use FastaReader;
use BowtieParser;
use DZUtil qw/overlap reverse_complement/;

my @range;
my $result = GetOptions (
    "range|r=i{2}" => \@range,
    "seqid|s=s" => \(my $seqid),
    "bisulfite|bs" => \(my $bs),
    "reference|f=s" => \(my $reference),
);
pod2usage(-verbose => 2, -noperldoc => 1) if ! @range || !$seqid || !$reference;

my $p = BowtieParser->new(file => \*ARGV);
my $fr = FastaReader->new(file => $reference, slurp => 1);

# PARSED_MATCHES is [ABS_COORD, ORIGINAL, NEW]

while (defined(my $bt = $p->next(1))){
    my ($readid, $strand, $chr, $read_pos, $read, $quality, $foobar, $mm_string, $mismatches) = @$bt;
    my $read_end = $read_pos + length($read) - 1;
    my $read_length = length($read);

    my $fix_rc = ($bs && $chr =~ s/^RC_//);
    if ($fix_rc){
        $strand = '-';
        ($read_pos, $read_end) = (
            $fr->reverse2forward($chr, $read_end),
            $fr->reverse2forward($chr, $read_pos),
        );
    }

    next if uc $chr ne uc $seqid;
    next if overlap([$read_pos, $read_end], \@range) == 0;

    @$mismatches = sort{
        $a->[0] <=> $b->[0]
    }
    map {
        my ($pos, $original, $new) = @$_;
        if ($fix_rc){
            $pos = $fr->reverse2forward($chr, $pos);
            $original = reverse_complement $original;
            $new = reverse_complement $new;
        }

        $range[0] <= $pos && $pos <= $range[1] ? ([$read_length - ($pos - $read_pos) - 1, $original, $new]) : ();
    } @$mismatches;

    if (@$mismatches){
        say join "\t", $readid, $strand, $chr, $read_pos, $read, $quality, $foobar, join ",", map { "$_->[0]:$_->[1]>$_->[2]" } @$mismatches; 
    }
}

=head1 NAME

bowtie-extract-range.pl - extract bowtie reads in a region of interest

=head1 SYNOPSIS

Usage examples:

Extract all reads overlapping 10000-20000 of chr4

 bowtie-extract-range.pl -r 10000 20000 -s chr4 -f TAIR_reference.fas input.bowtie > output.bowtie

extract a range from bs-bowtie (c2t reads vs c2t genome) output, converting
RC_* sequence coordinates appropriately:

 bowtie-extract-range.pl -bs -r 10000 20000 -s chr4 -f TAIR_reference.fas input.bowtie > output.bowtie

=cut
