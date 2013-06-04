#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Pod::Usage;
use Getopt::Long;
use IO::File;
use FindBin;
use lib "$FindBin::Bin/lib";
use BowtieParser;
use BigArray;
use FastaReader;
use Eland::Parser;
use GFF::Parser;
use Sam::Parser;
END {close STDOUT}
$| = 1;

my $result = GetOptions (
    "reference|r=s"   => \my $reference,
    "window-size|w=i" => \(my $window_size = 50),
    "no-skip|k"       => \(my $noskip),
    # "base|b=i"        => \(my $base = 1),
    "format|f=s"      => \(my $format),
    "strand|s"        => \(my $do_strand),
    "verbose|v"       => \(my $verbose),
    "check-rc|rc"     => \(my $check_rc),
    "output|o=s" => \(my $output = '-'),
    "first-base-only|1" => \(my $first_base_only),
);

pod2usage(-verbose => 2, -noperldoc => 1) 
if (!$result || !$reference || ! $format 
    || $format !~ /^(?:gff|eland|bowtie|sam|g|e|b|s)$/ 
    || ($format !~ /^(?:sam|s)$/ && $first_base_only));

my $fasta_reader = FastaReader->new(file => $reference, slurp => 0);

my %counters = map { 
    uc($_) => {
        (
            $do_strand ?  (
                '+' => BigArray->new(base => 1, size => $fasta_reader->get_length($_)),
                '-' => BigArray->new(base => 1, size => $fasta_reader->get_length($_))
            ) : (
                '.' => BigArray->new(base => 1, size => $fasta_reader->get_length($_)),
            )
        ),
    }
} $fasta_reader->sequence_list();

my %touched; # record touched sequences in output only

{
    my $c = 0;
    sub counter { 
        say STDERR $c if ($verbose && ++$c % 50000 == 0); 
    }
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
	    $touched{uc $chr}++;
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
        $touched{uc $gff->sequence}++;
    }
}
elsif ($format eq 'b' || $format eq 'bowtie'){
    my $bowtie_reader = BowtieParser->new(file => \*ARGV);
    while (defined(my $bowtie = $bowtie_reader->next())){
        my (undef, $strand, $chr, $pos, $read) = @$bowtie;
        my $end = $pos + length($read) - 1;
        $strand //= '+';

        if ($check_rc and $chr =~ s/^RC_//){
            ($end, $pos) = $fasta_reader->range_reverse2forward($chr, $pos, $end);
            $strand = $strand eq '+' ? '-' : '+';
        }
        if (! exists $counters{uc $chr}){
            if ($chr =~ /^RC_/){
                die "can't find seqid $chr in $reference, did you need to do -rc?";
            }
            die "can't find seqid $chr in $reference";
        }

        if ($do_strand){
            $counters{uc $chr}{$strand}->increment_range($pos, $end);
        }
        else{
            $counters{uc $chr}{'.'}->increment_range($pos, $end);
        }
        counter();
        $touched{uc $chr}++;
    }
}
elsif ($format eq 's' || $format eq 'sam'){
    my $sam_reader = Sam::Parser->new(file => \*ARGV, convert_rc => $check_rc, skip_unmapped => 1);
    while (defined(my $sam = $sam_reader->next())){
        my $strand = $sam->is_reverse() ? '-' : '+';
        my $first_base = $sam->is_reverse() ? $sam->rightmost : $sam->leftmost;

        if ($do_strand && $first_base_only){
            $counters{uc $sam->seqid}{$strand}->push_increment($first_base);
        }
        elsif ($do_strand && ! $first_base_only){
            for my $aln ($sam->matched_chunks) {
                $counters{uc $sam->seqid}{$strand}->increment_range(@$aln);
            }
            # $counters{uc $sam->seqid}{$strand}->increment_range($sam->leftmost, $sam->rightmost);
        }
        elsif (! $do_strand && $first_base_only){
            $counters{uc $sam->seqid}{'.'}->push_increment($first_base);
        }
        else { # ! $do_strand && ! $first_base_only
            for my $aln ($sam->matched_chunks) {
                $counters{uc $sam->seqid}{'.'}->increment_range(@$aln);
            }
            # $counters{uc $sam->seqid}{'.'}->increment_range($sam->leftmost, $sam->rightmost);
        }
        counter();
        $touched{uc $sam->seqid}++;
    }
}

say STDERR "Done counting, now outputting to $output" if $verbose;
say STDERR Dumper(\%touched) if $verbose;

my @strands = $do_strand ? qw/+ -/ : qw/./;

my $output_fh = $output eq '-' ? *STDOUT : IO::File->new($output, 'w');

# warn Dumper [sort $fasta_reader->sequence_list()];

for my $seq (sort $fasta_reader->sequence_list()) {
	unless (exists($touched{uc $seq}) || $noskip){ 
        # warn "skipping $seq" ; 
        next;
    }
    my $numreads = $touched{uc $seq};
    say STDERR "outputting $seq ($numreads reads)" if $verbose;

    my $start = 1;
    my $max = $fasta_reader->get_length($seq);
    while ($start <= $max){
        if ($window_size == 1){
            for my $s (@strands) {
                my $value  = $counters{uc $seq}{$s}->get_pdl()->at($start - 1); 
                if ($value > 0 || $noskip){
                    $output_fh->print(join "\t", $seq, qw/. ./, $start, $start, $value, $s, qw/. ./);
                    $output_fh->print("\n");
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
                    $output_fh->print(join "\t", $seq, qw/. ./, $start, $end, $count , $s, qw/. ./); 
                    $output_fh->print("\n");
                }
            }
            $start += $window_size;
        }
    }
}

close $output_fh if $output ne '-';

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

Format of alignment file. Can be "gff", "g", "eland", "e", "bowtie", "b", "sam", "s"  

=item --window-size <window_size> | -w <window_size>

Default 50.

=item --no-skip | -k 

Print window even if nothing maps to it.

=item --strand | -s 

Preserve strand information.  Default off.

=item --verbose | -v 

=item --first-base-only | -1 

Only supported for SAM file format currently.

=back

=cut

