#!/usr/bin/env perl
use Data::Dumper;
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib";
use FastaReader;
use Test::More qw(no_plan);
use feature 'say';

sub rc{
    my $seq = shift;
    $seq =~ tr/acgtACGT/tgcaTGCA/;
    return reverse $seq;
}

for my $slurp (0, 1) {
    my $f = FastaReader->new(file => 't/test.fasta', ht => sub { s/>(\w+)/$1/; return $_ }, slurp => $slurp );
    my %lengths = $f->sequence_lengths();

    #######################################################################
    # Basics
    
    is_deeply($f->length, { CHR1 => 711, CHR2 => 1185, CHR3 => 711, CHR4 => 474, CHR5 => 711, CHRC => 711, CHRM => 1343,}, "seq size with internal length hash");
    is_deeply({$f->sequence_lengths}, { chr1 => 711, chr2 => 1185, chr3 => 711, chr4 => 474, chr5 => 711, chrc => 711, chrm => 1343,}, "seq size with sequence_length");
    is_deeply([sort $f->sequence_list], [sort qw/chr1 chr2 chr3 chr4 chr5 chrc chrm/], "sequence_list");

    for ($f->sequence_list){
        ok($f->has_sequence($_), "has_sequence");
    }
    for (1 .. 10){
        ok(! $f->has_sequence('khr' . int(rand 100)), "! has_sequence");
    }


    #######################################################################
    # Get whole
    
    while (my ($seq,$len) = each %lengths) {
        is(length($f->get($seq,undef, undef)), $len, "whole seq length for $seq");
        if ($slurp){
            is($f->get($seq, 1, $lengths{$seq}, coord => 'f' , base => 1) , $f->get($seq) , "get($seq) whole implied vs get($seq) whole explicit");
        }
    }

    #######################################################################
    # Subseq
    
    is($f->get('CHR1' , 1  , 10 , coord => 'f' , base => 1) , 'CCCTAAACCC' , "base 1 forward 1-10 (slurp=$slurp)");
    is($f->get('CHr1' , 11 , 20 , coord => 'f' , base => 1) , 'TAAACCCTAA' , "base 1 forward 11-20 (slurp=$slurp)");
    is($f->get('ChR1' , 0  , 9  , coord => 'f' , base => 0) , 'CCCTAAACCC' , "base 0 forward 0-9 (slurp=$slurp)");
    is($f->get('chr1' , 10 , 19 , coord => 'f' , base => 0) , 'TAAACCCTAA' , "base 0 forward 10-19 (slurp=$slurp)");

    is($f->get('CHR1' , 1  , 10 , coord => 'f' , rc => 1, base => 1) , rc('CCCTAAACCC') , "base 1 forward rc 1-10 (slurp=$slurp)");
    is($f->get('CHr1' , 11 , 20 , coord => 'f' , rc => 1, base => 1) , rc('TAAACCCTAA') , "base 1 forward rc 11-20 (slurp=$slurp)");
    is($f->get('ChR1' , 0  , 9  , coord => 'f' , rc => 1, base => 0) , rc('CCCTAAACCC') , "base 0 forward rc 0-9 (slurp=$slurp)");
    is($f->get('chr1' , 10 , 19 , coord => 'f' , rc => 1, base => 0) , rc('TAAACCCTAA') , "base 0 forward rc 10-19 (slurp=$slurp)");

    is($f->get('ChR1' , 1  , 10 , coord => 'r' , base => 1) , rc('TTTTAGATGT') , "base 1 reverse rc 1-10 (slurp=$slurp)");
    is($f->get('CHR1' , 11 , 20 , coord => 'r' , base => 1) , rc('AAAAAAGTAT') , "base 1 reverse rc 11-20 (slurp=$slurp)");
    is($f->get('CHr1' , 0  , 9  , coord => 'r' , base => 0) , rc('TTTTAGATGT') , "base 0 reverse rc 1-10 (slurp=$slurp)");
    is($f->get('cHR1' , 10 , 19 , coord => 'r' , base => 0) , rc('AAAAAAGTAT') , "base 0 reverse rc 11-20 (slurp=$slurp)");

    is($f->get('ChR1' , 1  , 10 , coord => 'r' , rc => 0, base => 1) , 'TTTTAGATGT', "base 1 reverse 1-10 (slurp=$slurp)");
    is($f->get('CHR1' , 11 , 20 , coord => 'r' , rc => 0, base => 1) , 'AAAAAAGTAT', "base 1 reverse 11-20 (slurp=$slurp)");
    is($f->get('CHr1' , 0  , 9  , coord => 'r' , rc => 0, base => 0) , 'TTTTAGATGT', "base 0 reverse 1-10 (slurp=$slurp)");
    is($f->get('cHR1' , 10 , 19 , coord => 'r' , rc => 0, base => 0) , 'AAAAAAGTAT', "base 0 reverse 11-20 (slurp=$slurp)");

    my $chr1len = $lengths{chr1};
    is($f->get('CHR1' , 1  , 1 , coord => 'f' , base => 1) , 'C' , "base 1 forward 1-1 (slurp=$slurp)");
    is($f->get('CHR1' , $chr1len, $chr1len, coord => 'f' , base => 1) , 'T' , "base 1 forward last-last (slurp=$slurp)");

    is($f->get('chr1' , $chr1len - 10, $chr1len-1, coord => 'f' , base => 0) , 'TTTTAGATGT' , "base 0 forward last ten (slurp=$slurp)");
    is($f->get('cHr1' , $chr1len - 9, $chr1len, coord => 'f' , base => 1) , 'TTTTAGATGT' , "base 1 forward last ten (slurp=$slurp)");

}



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

