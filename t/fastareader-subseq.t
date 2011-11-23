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
#for my $slurp (0) {
    my $f = FastaReader->new(file => 't/data/test.fasta', ht => sub { s/>(\w+)/$1/; return $_ }, slurp => $slurp );
    my %lengths = $f->sequence_lengths();
    my $chr1len = $f->get_length("chr1");
    my $chr4len = $f->get_length("chr4");
    my $chrmlen = $f->get_length("chrm");

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
    
    say Dumper \%lengths;
    while (my ($seq,$len) = each %lengths) {
        is(length($f->get($seq,undef, undef)), $len, "whole seq length for $seq");
        if ($slurp){
            is($f->get($seq, 1, $lengths{$seq}, coord => 'f' , base => 1) , $f->get($seq) , "get($seq) whole implied vs get($seq) whole explicit");
        }
    }
    is(substr($f->get('chr1', undef, undef, coord => 'f', base => 1,rc => 1),0, 10), rc('TTTTAGATGT'), "get whole implicit rc");

    #######################################################################
    # Subseq
    
    is($f->get('CHR1' , 1  , 10 , coord => 'f' , base => 1) , 'CCCTAAACCC' , "base 1 forward 1-10 (slurp=$slurp)");
    is($f->get('CHr1' , 11 , 20 , coord => 'f' , base => 1) , 'TAAACCCTAA' , "base 1 forward 11-20 (slurp=$slurp)");
    is($f->get('ChR1' , 0  , 9  , coord => 'f' , base => 0) , 'CCCTAAACCC' , "base 0 forward 0-9 (slurp=$slurp)");
    is($f->get('chr1' , 10 , 19 , coord => 'f' , base => 0) , 'TAAACCCTAA' , "base 0 forward 10-19 (slurp=$slurp)");

    for my $base (0,1){
        is($f->get('CHRm' , 1028-$base, 1037 -$base, coord => 'f' , base => 1-$base) , 'ACACCGTTTT' , "base $base forward beyond 1000 (slurp=$slurp)");
        is($f->get('CHrm' , 1028-$base, 1037 -$base, coord => 'f' , rc => 1,base => 1-$base) , rc('ACACCGTTTT') , "base $base forward beyond 1000 rc (slurp=$slurp)");
        is($f->get('CHrm' , $chrmlen-9-$base, $chrmlen-$base, coord => 'r' , rc => 0,base => 1-$base) , 'GGATCCGTTC', "base $base reverse beyond 1000 (slurp=$slurp)");
        is($f->get('CHrm' , $chrmlen-9-$base, $chrmlen-$base, coord => 'r' , rc => 1,base => 1-$base) , rc('GGATCCGTTC') , "base $base reverse beyond 1000 rc (slurp=$slurp)");
    }
    
    is($f->get('CHrm' , 948, 1026, coord => 'f' , rc => 0,base => 0) , 
        'AAGTGAGCTCACTAGCTGCTTGTTGTCTTGTAGAGTAGAAGACTTATAGATTAAAATTCTCCAACATATAGATGTCCTT', 
        "base 0 across increment(slurp=$slurp)");

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

    is($f->get('CHR1' , 1  , 1 , coord => 'f' , base => 1) , 'C' , "base 1 forward 1-1 (slurp=$slurp)");
    is($f->get('CHR1' , $chr1len, $chr1len, coord => 'f' , base => 1) , 'T' , "base 1 forward last-last (slurp=$slurp)");

    is($f->get('chr1' , $chr1len - 10, $chr1len-1, coord => 'f' , base => 0) , 'TTTTAGATGT' , "base 0 forward last ten (slurp=$slurp)");
    is($f->get('cHr1' , $chr1len - 9, $chr1len, coord => 'f' , base => 1) , 'TTTTAGATGT' , "base 1 forward last ten (slurp=$slurp)");

    is($f->get('ChR1' , undef  , 10 , coord => 'r' , rc => 0, base => 1) , 'TTTTAGATGT', "base 1 reverse 1-10 first ten with undef start(slurp=$slurp)");
    is($f->get('ChR1' , undef  , 9  , coord => 'f' , rc => 1, base => 0) , rc('CCCTAAACCC') , "base 0 forward rc 0-9 with undef start(slurp=$slurp)");
    is($f->get('cHr1' , $chr1len - 9, undef, coord => 'f' , base => 1) , 'TTTTAGATGT' , "base 1 forward last ten with undef end(slurp=$slurp)");
    is($f->get('cHr1' , $chr1len - 4, undef, coord => 'f' , base => 1) , 'GATGT' , "base 1 forward last ten with undef end(slurp=$slurp)");

    #######################################################################
    # lenient

    is($f->get('CHR1' , 0  , 10 , coord => 'f' , base => 1, lenient => 1) , 'CCCTAAACCC' , "lenient base 1 forward 1-10 (slurp=$slurp)");
    is($f->get('ChR1' , -10  , 10 , coord => 'r' , base => 1, lenient => 1) , rc('TTTTAGATGT') , "lenient base 1 reverse rc 1-10 (slurp=$slurp)");
    is($f->get('cHr1' , $chr1len - 9, $chr1len + 1, coord => 'f' , base => 1,lenient => 1) , 'TTTTAGATGT' , "lenient base 1 forward last ten (slurp=$slurp)");

    #######################################################################
    # death
    #is($f->get('CHR1' , 0  , 0 , coord => 'f' , base => 1, lenient => 1) , 'CCCTAAACCC' , "lenient base 1 forward 1-10 (slurp=$slurp)");
    
    for my $base (0,1){
        is($f->get_context('chr1', 2+$base, base => $base, rc => 0), 'CHH', "context 1 base $base");
        is($f->get_context('chr1', 7+$base, base => $base, rc => 0), 'CHH', "context 2 base $base");
        is($f->get_context('chr2', 474+$base, base => $base, rc => 0), 'CG', "context 3 base $base");
        is($f->get_context('chr1', 160+$base, base => $base, rc => 0), 'CG', "context 4 base $base");

        is($f->get_context_raw('chr1', 2+$base, base => $base, rc => 0), 'CTA', "context_raw 1 base $base");
        is($f->get_context_raw('chr1', 7+$base, base => $base, rc => 0), 'CCC', "context_raw 2 base $base");
        is($f->get_context_raw('chr1', 160+$base, base => $base, rc => 0), 'CGT', "context_raw 4 base $base");

        # out of bounds should return CHH
        is($f->get_context('chr4', 1, rc => 1),'CHH', "context edge case 1");
        is($f->get_context('chr4', 2, rc => 1),'CHH', "context edge case 2");
        is($f->get_context('chr1', $chr1len, rc => 0),'CHH', "context edge case 3");
        is($f->get_context('chr4', $chr4len, rc => 0),'CHH', "context edge case 4");

        is($f->get_context_raw('chr4', 1, rc => 1),'C', "context_raw edge case 1");
    }
    
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
>chr2
CTGAAATAAGAGGAGTTCTTCCAGAATCGCGTCATGGCCTCGAAGACAGAGTCAAGCTTGGTGAAAAGCAGATCGTTCG
AGGAGCTCTGACCACGGCCGGCTACGAGATAGACCACCACGAACCGTCTTGGTCACACTGTTGACAAATGCATGGACGG
TTGGAGGTTGTGATTCCACCACATCAGTTACAATCGGAGGCTGATCCGGACTTGGCTGTTGCGGTGGTGGAAATCGGAA
GCACAGGGGGATTTGTGGACACCATTCGATTTTGTTCTTTAGGGATTTCAAGATTGGAAAATTTTAGGAGACGATTGAT
TTGGAAGTCGGTGAAGTTTCGGTAGTGTTATGTATTTATTTTTTTTGTCATAAGCGGTAGCGTTATGTATTTGTAAGAT
CACATTAAGCGTCTTTCAACTTTGATTATACGCATTGCTTCATCTGATATCTAATCCATATTTTGATCAAAGTGACATT
TGATAAAATAAAATAAAACTAACCTATTAATATTATAACTGAGGTTATTTTTTTAGAAAAAAGATAGCTAGCTACACGA
ATTATACGAGAATACATGGTCACATGCCAAGTATATTATTGTATCCTTAACTACACGAAATTTATGTATTGCATGGCCA
CATATATGAACTATACTACGTCATGTGGCTAACATCATAGACATATAATTCCGGTATATCGACTGGTTGACTCTGGCTT
TGACTAACATTGACCGGCGTTGACCAACAAAAAATTTCAGAAAAAAACTTTAAAATAGTTTTTAATATTAAAAAATAAG
AAACTGTTTTTAGGTTTTGTATAAGAAAAAAAATTCTATTTTCAGCATAGATAGATTTGTATTTTTATCCATTAAATTT
TATAATAATTAGTAAAAAGTCATTATTAAATTTTAAAATATTAAAATGGAAAAATATTATTTGAACAATTAAGTAGGAA
ATTGCTTAAATTTATATTGTCAATTAAAAAACCTTAATTTCTATACTATTTTTATTATTATTTTGTATTTCTCAACCAC
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>chr3
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNCCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACC
CTAAATCCATAAATCCCTAAAACCATAATCCTAAATCCCTTAATTCCTAAATCCCTAATACTTAGACCCTAATCTTTAG
TTCCTAGACCCTAATCTTTAGTTCCTAGACCCTAAATCCATAATCCTTAATTCCTAAATTCCTAAATCCCTAATACTAA
ATCTCTAAATCCCTAGCAATTTTCAAGTTTTGCTTGATTGTTGTAGGATGGTCCTTTCTCTTGTTTCTTCTCTGTGTTG
TTGAGATTAGTTTGTTTAGGTTTGATAGCGTTGATTTTGGCCTGCGTTTGGTGACTCATATGGTTTGATTGGAGTTTGT
TTCTGGGTTTTATGGTTTTGGTTGAAGCGACATTTTTTTGTGGAATATGGTTTTTGCAAAATATTTTGTTCCGGATGAG
TAATATCTACGGTGCTGCTGTGAGAATTATGCTATTGTTTTGCAGGTCCTGTTCTTAATCTTTCATCGCTTTTGTGCTT
ATTGTCTCCTTGTCGTTTATGTTGAGTGGTGTTTGGGCTTTAGTCAGGTTCAGTCCAACCGTTGCTCTTCCTCCCGTTC
>chr4
GGAATGGTCGACGAAAATCTATCGGGTTCGAGGATTCGTCGACCAGGTGTTGGAATCGTCGACCGAGTCTGAGAATTCG
TAGACCAGGACGGCGGAATCCTCGACAATGACGAGGTATGGTCGAGGAAAATCTATCGGGTTCGAGGATTCGTCTACCA
GGTGATGGAATCCTCGACCAGGACAAAGAATTCGTCGACCAGGGGTGGAATTGTTGTTATTCCGATCATGAGAGCGGAT
ATCAGTACAGATCCGACGCTGGTGAAAAAGATCACGGCGATCGTGGATAGTATCAAGCCACCGAGAGTCTCGTATTCGG
AGAAAGATCGGCCGATGAGGATAAGAGGTCGATCGGATGGACGGAAGAGGTAGAGGAAGAGCCATGAAGCGGCGAGGCA
TAGGAGGAGGATGAGCGAGAATGGGTGGGCGGGAAGAGAGAAACTGATGATCAGAGCGATGATGCAGACGTAATTCACC
>chr5
TATACCATGTACCCTCAACCTTAAAACCCTAAAACCTATACTATAAATCTTTAAAACCTATACTCTAAACCATAGGGTT
TGTGAGTTTGCATAAAGTGTCACGTATAAGTGTTTCTAACATGTGAGTTTGCATAAGAGTCTCGACTATGTGTTTGTTC
AAAAGTGACGTAAGTGTTTAGACTAGAGCCGGCCGTGAGCACAAGCGGGCCAAGCCCATGCTTGCGGCAGATTATCCTA
TAATTATGTTTTGCGGCTTTACAATTTTGAATTTGTTTTGTTGTTTAGTGGTCGAGTCAGGGATGAGTTTGTTCTCAAA
CATCTCAAATTTCTAATCTTCCAAGTGTTAGGTTCACTCTACCTTATTTTTTTTTTTTTTTTTTATTTTAAATTTCTAT
TAGAAAGAACTTGAACCTTATTTCAAGTGAACATTAAAAAGAATTGTAACCCTACCTTATTTCTATTATTATGAATATG
TGTAAAGTAAATTTTTTTACAATAAGCATGTATAAAGTAACATATTTTTGTTCTTTTCTTATAGCCTTGTAAATCATAG
GGACGGATCTAGTTTAGAATAAAGTGTTGCATAAGTTCTTCGACTAACATGTGAGTCCATTTTTGTGAGTTTTGCATAA
CTTTTGCGCCTTCGACAATTTGAATTTTATGCATAAGTGCTTTGACTAATTTATATGTTTTGTCTAAGTGATTCGACAT
>chrc
ATGGGCGAACGACGGGAATTGAACCCGCGATGGTGAATTCACAATCCACTGCCTTAATCCACTTGGCTACATCCGCCCC
TACGCTACTATCTATTCTTTTTTGTATTGTCTAAAAAAAAAAAAAAATACAAATTTCAATAAAAAATAAAAAAAGGTAG
CAAATTCCACCTTATTTTTTTTCTAATAAAAAATATATAGTAATTTTTTATTATTTATTATTATTATTTATTATTAATA
TAATAAATAAAGTAAAATATGATACTCTATAAAAATTTGCTCATTTTTATAGAAAAAAACGAGTAATATAAGCCCTCTT
TCTTATTTAAAGAAGGCTTATATTGCTCGTTTTTTACTAAACTAGATCTAGACTAACACTAACGAATTATCCATTTGTA
GATGGAGCCTCAACAGCAGCTAGGTCTAGAGGGAAGTTGTGAGCATTACGTTCATGCATAACTTCCATACCAAGGTTAG
CACGGTTAATAATATCAGCCCAAGTATTAATAACACGTCCTTGACTATCAACTACTGATTGGTTGAAATTGAAACCATT
TAGGTTGAAAGCCATAGTACTAATACCTAAAGCAGTAAACCAAATACCTACTACCGGCCAAGCCGCTAAGAAGAAATGT
AAAGAACGAGAATTGTTGAAACTAGCATATTGGAAAATCAATCGGCCAAAATAACCGTGAGCAGCTACAATGTTGTAAG
>chrm
GGATCCGTTCGAAACAGGTTAGCCTACTATAATATAAGGATTGGATTCTAATAAGTTCGAAACAGGTTAGCCTTAGCCT 0
ACTATAGGATTAGATCTTTCTTATCAACCTACTAACTTCTTCCTTGTTGGGATGAGAAACCCTTTTGCAACCAAGCGTG 79
CTTTGAGTTTGTCAAGGGACCCATCTGCATTCAGTTTCACTCTGAAAACCCATTTACAACCGAGAAGATTCATGTCAGG 158
TGATGCGGGAACTAAGTCCCAAGTGTGATTCTGTGTTAATGCCGACATCTCTTCTTGCATAGCTTGTCTCCATCCTGGG 237
AGGCAGACGTAATGGTTTTTGGTTCAGAGGGAGTGTATTTTTGTGTAAACAGGTTGTAACGAGGATTAGGCTTGCGAAT 316
ACCATCCTTTGCCCGAGTGATCATATGATGTCTATTAGGTGAAAGTAGCTCAGGAGCAGCTGTCCCAACATCAAAAAAG 395
GTACCGCTGTCGCCAATAGGAACAGGATCTGAGCCTGCCGTACGCACAGGACAGTCTCTTTCTGATGTGGTAGCAGTTC 474
CAGGAGCAATCGCAGAGACAATTGAAGGATCTGCAGATTCGCAATCAGAGCCTGTGAAACCGGGAAGATGTCGAGATAG 553
CAGGAGAAGAATTAGCAGTATGTGAGTGCGGAAGCAGAGTGGAGGATCGAGGGACCTGTGAAGGTGTGATAAAGGAACT 632
GTCTGTTGAAACAGAAGGGACATATGAAGGTTTAGAACCAAGCTGCCAAGCACGAAGCAGAGAACTTGGATACGGTGGG 711
ATGAGATGGTGATAACAGTCCTGAAAACGAAAGACTAGAGCTCTTGCCCAACTGACAGGCTTCACAAACAGAAAGATCC 790
CTTTTATTAATAACTATAGCTTTACTAGTCTTGAGTTGTTGGAGAACTTGAGGATTTGGATGCCCAAGCCTATAATGCC 869
AAGTGAGCTCACTAGCTGCTTGTTGTCTTGTAGAGTAGAAGACTTATAGATTAAAATTCTCCAACATATAGATGTCCTT 948
ACACCGTTTTCCTTTGCTCAGCAGGCTCCGTGTTTGCTTGTCCTTTATGCATACTTCGTTAGCATCAAACTCAAAGAAG 1027
CATGGATAAATCATCACAAAGTTTGGACACAGAGAGAAGAGACTTAGTTATGAAAGGACAAACAAGAACTTCATTTAAT
GGTAAACTACCTGATGAGCTTGGTAGATTGGTAGATCCAACATGAGTGATGGGCAAGAAGGCTCCATCCCCAACCATCA
CACTGTCGGATCCCACATAAGGTTGTGACTGTTGGAGGTGATGAGCAGAATTGGTAATGTGAGAGGAGGCACCAGAATC
