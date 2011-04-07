#!/usr/bin/env perl
use Data::Dumper;
use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/../lib";
use Fasta;

use Test::More qw(no_plan);
use feature 'say';

my $seq = slurp_fasta('t/test.fasta');

is(fasta_subseq($seq , 'CHR1' , 1  , 10 , reverse => 0 , base => 1) , 'CCCTAAACCC' , "base 1 forward 1-10");
is(fasta_subseq($seq , 'CHr1' , 11 , 20 , reverse => 0 , base => 1) , 'TAAACCCTAA' , "base 1 forward 11-20");
is(fasta_subseq($seq , 'ChR1' , 0  , 9  , reverse => 0 , base => 0) , 'CCCTAAACCC' , "base 0 forward 0-9");
is(fasta_subseq($seq , 'chr1' , 10 , 19 , reverse => 0 , base => 0) , 'TAAACCCTAA' , "base 0 forward 10-19");

sub rc{
    my $seq = shift;
    $seq =~ tr/acgtACGT/tgcaTGCA/;
    return reverse $seq;
}
say length $seq->{'CHR1'};
is(fasta_subseq($seq , 'ChR1' , 1  , 10 , reverse => 1 , base => 1) , rc('TTTTAGATGT') , "base 1 reverse 1-10");
is(fasta_subseq($seq , 'CHR1' , 11 , 20 , reverse => 1 , base => 1) , rc('AAAAAAGTAT') , "base 1 reverse 11-20");
is(fasta_subseq($seq , 'CHr1' , 0  , 9  , reverse => 1 , base => 0) , rc('TTTTAGATGT') , "base 0 reverse 1-10");
is(fasta_subseq($seq , 'cHR1' , 10 , 19 , reverse => 1 , base => 0) , rc('AAAAAAGTAT') , "base 0 reverse 11-20");

__DATA__
# test.fasta:
>chr1
CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCAT
GAATCCCTAAATACCTAATTCCCTAAACCCGAAACCGGTTTCTCTGGTTGAAAATCATTGTGTATATAATGATAATTTT
ATCGTTTTTATGTAATTGCTTATTGTTGTGTGTAGATTTTTTAAAAATATCATTTGAGGTCAATACAAATCCTATTTCT
TGTGGTTTTCTTTCCTTCACTTAGCTATGGATGGTTTATCTTCATTTGTTATATTGGATACAAGCTTTGCTACGATCTA
CATTTGGGAATGTGAGTCTCTTATTGTAACCTTAGGGTTGGTTTATCTCAAGAATCTTATTAATTGTTTGGACTGTTTA
TGTTTGGACATTTATTGTCATTCTTACTCCTTTGTGGAAATGTTTGTTCTATCAATTTATCTTTTGTGGGAAAATTATT
TAGTTGTAGGGATGAAGTCTTTCTTCGTTGTTGTTACGCTTGTCATCTCATCTCTCAATGATATGGGATGGTCCTTTAG
CATTTATTCTGAAGTTCTTCTGCTTGATGATTTTATCCTTAGCCAAAAGGATTGGTGGTTTGAAGACACATCATATCAA
AAAAGCTATCGCCTCGACGATGCTCTATTTCTATCCTTGTAGCACACATTTTGGCACTCAAAAAAGTATTTTTAGATGT

