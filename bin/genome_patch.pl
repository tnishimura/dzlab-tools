#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;

use FindBin;
use lib "$FindBin::Bin/../lib";
use Fasta;

use Pod::Usage;
use Getopt::Long;

my $help;
my $reference;
my $result = GetOptions (
    "reference|r=s" => \$reference,
    "allow-mismatch|m" => \(my $allow_mismatch),
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
        if ($grab ne $original){
            ++$counter;
            if ($allow_mismatch){
                say STDERR "BEFORE: $grab doesn't match $original, but converting anways"; 
                say STDERR "line was: $line";
            }
            else{
                say STDERR "BEFORE: $grab doesn't match $original, skipping SNP"; 
                say STDERR "line was: $line";
                next TOP;
            }
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

By default, if the "from" base in column 9 does not match what is found in the
genome, an error will be printed and that SNP will be skipped. With
--allow-mismatch/-m, the base will be converted to the "to" base even if there
is an error.

 genome_patch.pl -m -r reference.fas snplist.gff > new-reference.fas

=cut
