#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;

use FindBin;
use lib "$FindBin::Bin/lib";
use Fasta;

use Pod::Usage;
use Getopt::Long;

my $help;
my $reference;
my $result = GetOptions (
    "reference|r=s" => \$reference,
);
if (! $result || ! $reference || (! @ARGV && -t STDIN)){
    pod2usage(-verbose => 2, -noperldoc => 1);
}

my $genome = slurp_fasta($reference);

my $counter = 0;
TOP:
while (defined(my $line = <ARGV>)){
    $line =~ tr/\r\n//d;
    my ($seq, $coord, $strand,$attr) = (split /\t/, $line)[0,3,6,8];
    
    next unless $seq && $coord && $attr;
    my ($original, $converted) = split />/,$attr;
    die "$line fubar" unless $original && $converted;
    if ($strand eq '-'){
        $original =~ s/acgtACGT/tgcaTGCA/;
        $converted =~ s/acgtACGT/tgcaTGCA/;
    }

    $seq = uc $seq;

    # check before
    BEFORE:
    {
        my $grab = substr $genome->{$seq}, $coord-1,1;
        last BEFORE if $grab !~ /ACTG/ && $original !~ /ACTG/;
        if ($grab ne $original){
            say STDERR "BEFORE: $grab doesn't match $original at $line"; 
            ++$counter;
            next TOP;
        }
    }

    substr($genome->{$seq}, $coord-1,1) = $converted;

    # check after
    {
        my $grab = substr $genome->{$seq}, $coord-1,1;
        warn "AFTER: $grab doesn't match $original at $line" unless $grab eq $converted;
    }
}

for my $seq (sort keys %$genome) {
    say format_fasta($seq, $genome->{$seq});
}
=head1 genome_patch.pl 

 genome_patch.pl -r reference.fas snplist.gff > new-reference.fas

snplist.gff should be a gff file with the following columns:

 column 1: sequence id 
 column 4: coord
 column 7: strand (if '.', assumed to be + strand).
 column 9: SNP, written as "A>C" or "C>G".  Should be with respect to
           the strand (column 7). 

=cut
