#!/usr/bin/env perl
use warnings;
use strict;
use 5.010_000;
use Data::Dumper;
use autodie;

use Getopt::Long;
use Pod::Usage;

use FindBin;
use lib "$FindBin::Bin/lib";
use FastaReader;

my @excluded_sequences;

# Grabs and parses command line options
my $result = GetOptions (
    'kmer|k=i'           => \(my $kmer = 2),
    'exclude-seq|x=s{,}' => \@excluded_sequences,
    'alphabet|a=s'       => \(my $alphabet = 'ACGT'),
    'output|o=s'         => \(my $output),
);

pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || ! @ARGV);

if ($output) {
    open my $USER_OUT, '>', $output;
    select $USER_OUT;
}

my $reference = shift;
my %nucleotide_frequencies;

say STDERR "# Reading in $reference";
my $fr = FastaReader->new(file => $reference, slurp => 0);
say STDERR "# Done";

for my $k (1, $kmer) {
    
    say STDERR "#k:\t$k";

    my %frequencies  = ();
    my $total_length = 0;

    for my $chromosome ($fr->sequence_list()) {
        my $chrlen = $fr->get_length($chromosome);

        $total_length += $chrlen;
        
        say STDERR "#chromosome:\t$chromosome";
        say STDERR "#length:\t", $chrlen;

        my %frequency = %{ word_composition ($fr->get($chromosome, undef, undef), $k) };

        for my $word (sort keys %frequency) {

            next if $word =~ m/[^$alphabet]/i;

            $frequencies{$word} += $frequency{$word};

            say STDERR join ("\t",
                               $word,
                               $frequency{$word},
                               $frequency{$word} / $chrlen, 
                           );
        }
    }

    say STDERR "#all:";
    say STDERR "#length:$total_length";

    say join ("\t",
                '#word',
                '#count',
                '#size',
                '#observed',
                '#expected',
                '#obs/expect',
                '#independent',
            )
            if $k > 1;

    for my $word (sort keys %frequencies) {

        next if $word =~ m/[^$alphabet]/i;

        if ($k == 1) {
            $nucleotide_frequencies{$word} = $frequencies{$word} / $total_length;
            say STDERR join ("\t",
                        $word,
                        $nucleotide_frequencies{$word},
                    );
        }
        else {
            my $observed = $frequencies{$word} / $total_length;
            my $expected = 1;
            map {$expected *= $nucleotide_frequencies{$_}} (split //, $word);

            say join ("\t",
                        $word,
                        $frequencies{$word},
                        $total_length,
                        $observed,
                        $expected,
                        $observed / $expected,
                        map {$nucleotide_frequencies{$_} * $total_length} (split //, $word),
                    );
        }
    }
}

sub word_composition {
    my ($sequence, $k, $alphabet) = @_;

    my %frequency = ();

    for (0 .. length ($sequence) - $k) {
     
        my $word = substr $sequence, $_, $k;
        next if $word =~ m/[^ACGT]/i;
        $frequency{$word}++;

    }

    return \%frequency;
}

=head1 NAME

nucleotide-freqs.pl - calculate the nucleotide and nucleotide k-mer
composition. 

=head1 SYNOPSIS

Calculate 1-mer and 2-mer composition in reference.fasta:

 nucleotide-freqs.pl -o output.txt reference.fasta 
 nucleotide-freqs.pl reference.fasta > output.txt

Calculate 3-mers as well:

 nucleotide-freqs.pl reference.fasta 

 nucleotide-freqs.pl -o output.txt --kmer 3 reference.fasta 
 nucleotide-freqs.pl -o output.txt -k     3 reference.fasta 

=head1 OPTIONS

=over

=item  -k <count> | --kmer <count>

Calculate composition up to this long. Default 2.  This means that A, C, G, T,
AA, AC, AG, AT, ..., TG, TT composition are calculated.

=item  -a <alpha> | --alphabet <alpha>

Alphabet of genome.  Default ACGT.  Use, say, ACGU for RNA.

=item  -o <file> | --output <file>

Output file. Default to standard out.

=back

=cut
