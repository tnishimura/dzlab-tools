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
use BowtieParser;
use BigArray;
use FastaReader;
use Eland::Parser;
use GFF::Parser;
END {close STDOUT}
$| = 1;

my $result = GetOptions (
    "reference|r=s"   => \my $reference,
    "window-size|w=i" => \(my $window_size = 50),
    "no-skip|k"       => \(my $noskip),
    "base|b=i"        => \(my $base = 1),
    "format|f=s"      => \(my $format),
    "strand|s"        => \(my $do_strand),
    "verbose|v"       => \(my $verbose),
    "output|o=s" => \(my $output = '-'),
);

pod2usage(-verbose => 2, -noperldoc => 1) 
if (!$result || !$reference || ! $format || $format !~ /^(?:gff|eland|bowtie|g|e|b)$/);

my $fasta_reader = FastaReader->new(file => $reference, slurp => 0);

my %counters = map { 
    uc($_) => {
        (
            $do_strand ?  (
                '+' => BigArray->new(base => $base, size => $fasta_reader->get_length($_)),
                '-' => BigArray->new(base => $base, size => $fasta_reader->get_length($_))
            ) : (
                '.' => BigArray->new(base => $base, size => $fasta_reader->get_length($_)),
            )
        ),
    }
} $fasta_reader->sequence_list();

{
    my $c = 0;
    sub counter { say STDERR $c if ($verbose && ++$c % 50000 == 0); }
}

if ($format eq 'e' || $format eq 'eland'){
    my $eland_parser = Eland::Parser->new(file => \*ARGV, fastareader => $fasta_reader);
    while (defined(my $eland = $eland_parser->next())){
        my (undef, undef, $positions) = @$eland;
        for my $pos (@$positions) {
            #say STDERR Dumper $pos;
            my ($chr, undef, $is_reverse, $start, $end) = @$pos;
            if ($do_strand){
                $counters{uc $chr}{$is_reverse ? '-' : '+'}->increment_range($start, $end);
            }
            else{
                $counters{uc $chr}{'.'}->increment_range($start, $end);
            }
            counter();
        }
    }
}
elsif ($format eq 'g' || $format eq 'gff'){
    my $gff_parser = GFF::Parser->new(file => \*ARGV);
    while (defined(my $gff = $gff_parser->next())){
        if ($do_strand){
            my $strand = defined($gff->strand()) && $gff->strand() eq '-' ? '-' : '+';
            $counters{uc $gff->sequence}{$strand}->increment_range($gff->start(), $gff->end());
        }
        else{
            $counters{uc $gff->sequence}{'.'}->increment_range($gff->start(), $gff->end());
        }
        counter();
    }
}
elsif ($format eq 'b' || $format eq 'bowtie'){
    my $bowtie_reader = BowtieParser->new(file => \*ARGV);
    while (defined(my $bowtie = $bowtie_reader->next())){
        my (undef, $strand, $chr, $pos, $read) = @$bowtie;
        if ($do_strand){
            $strand = defined($strand) && $strand eq '-' ? '-' : '+';
            $counters{uc $chr}{$strand}->increment_range($pos, $pos + length($read) - 1);
        }
        else{
            $counters{uc $chr}{'.'}->increment_range($pos, $pos + length($read) - 1);
        }
        counter();
    }
}

my @strands = $do_strand ? qw/+ -/ : qw/./;

if ($output ne '-'){
    open my $output_fh, '>', $output;
    select $output_fh;
}

for my $seq (sort $fasta_reader->sequence_list()) {
    my $start = 1;
    my $max = $fasta_reader->get_length($seq);
    #say STDERR "$seq from $start to $max";
    while ($start <= $max){
        if ($window_size == 1){
            for my $s (@strands) {
                my $value  = $counters{uc $seq}{$s}->{pdl}->at($start - $base); 
                if ($value > 0 || $noskip){
                    say join "\t", $seq, qw/. ./, $start, $start, $value, $s, qw/. ./;
                }
            }
            ++$start;
        }
        else{
            my $end = $start + $window_size - 1;
            if ($end > $max){
                $end = $max;
            }
            for my $s (@strands) {
                my $pdl  = $counters{uc $seq}{$s}->get_range($start, $end);
                my $count  = $pdl->max();

                if ($count  > 0 || $noskip){ 
                    say join "\t", $seq, qw/. ./, $start, $end, $count , $s, qw/. ./; 
                }
            }
            $start += $window_size;
        }
    }
}

=head1 NAME

 window_alignment.pl - given an alignment file, return a gff file of windows
 scores showing how many reads overlap that window.  

=head1 SYNOPSIS

Usage examples:

 window_alignment.pl -r genome.fasta -w 1 -f bowtie alignment.bowtie
 window_alignment.pl -r genome.fasta -k -w 50 -f bowtie alignment.bowtie
 window_alignment.pl -r genome.fasta -k -w 50 -f gff bowtie_converted_to_gff.gff
 window_alignment.pl -r genome.fasta -k -w 50 -f eland bowtie_converted_to_eland.eland3

=head1 OPTIONS

=over

=item --reference <fasta> | -r <fasta>

Reference genome file.

*WARNING* This script will consume approximately (4 * size_of_reference_genome)
bytes of memory.  If --strand is used, it will use double that. 

=item --format <f> | -f <f>

Format of alignment file. Can be "gff", "g", "eland", "e", "bowtie", "b".  

=item --window-size <window_size> | -w <window_size>

Default 50.

=item --no-skip | -k 

Print window even if nothing maps to it.

=item --base <b> | -b <b>

The base of the coordinates in the file (ie, is the first coordinate of the
chromosomes 0 or 1?)  This should match the -B options passed to bowtie.
Default 1.

=item --strand | -s 

Preserve strand information.  Default off.

=item --verbose | -v 

=back

=cut

