#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More;
use Test::Exception;
use Scalar::Util qw/looks_like_number/;

use Sam::Alignment;
use GFF::Parser;
use FastaReader;
use IO::All;

use TestUtils;

sub check_basics{
    my $name = shift;
    my $is_mapped = shift;
    my $line = join "\t", @_;

    my $sam = Sam::Alignment->new($line);
    my $p = "$name (basic)"; # prefix
    if ($is_mapped){
        ok($sam->mapped, "$p: should've mapped");
        ok(looks_like_number($sam->leftmost) && $sam->leftmost > 0, "$p: positive position");
        ok($sam->cigar, "$p: has a cigar string");
        is(ref $sam->cigar, 'ARRAY', "$p: cigar looks valid");
        isnt($sam->seqid, '*', "$p: has a seqid");

        is(
            $sam->cigarlength(), 
            $sam->readlength(),
            "$p: cigar length is readlength"
        );

        ok($sam->original_mismatch_string, "$p: has mismatch string");
    }
    else{
        ok(! $sam->mapped, "$p: shouldn't have mapped");
        is($sam->original_cigar_string, '*', "$p: doesn't have a cigar string");
        is($sam->seqid, '*', "$p: doesn't have a seqid");
    }
    is(length($sam->readseq), length($sam->readqual), "$name (basic): readseq and readqual same length");

    return $sam;
}

#######################################################################
# basics


check_basics("01", 0, qw{
        HS2:242:D17PPACXX:6:1101:1070:2241	
        4	*	0	0	*	*	0	0	
        ATTGTAAATATTAGGTATTTTGTTTGGGGGGAAGGTTATATTTAGTATTGGATTAATGATTTGGGTGATATGTTTTAGGTTTTTTTTTTTATGTATTGTA	
        IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	
        YT:Z:UU
    });

check_basics("02", 1, qw{
    HS2:242:D17PPACXX:6:1101:1387:2206      
    0       RC_pBIN-pROK2   3903    2       39M61S  *       0       0 
    ATGGAGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGGTTTGGGGGGGGTGG    
    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
    AS:i:54 XS:i:50 XN:i:0  XM:i:3  XO:i:0  XG:i:0  NM:i:3  MD:Z:4T4G8G20   YT:Z:UU
    });

#######################################################################
# forward snps

{
    my $sam = Sam::Alignment->new(join("\t", qw{
    HS2:242:D17PPACXX:6:1101:1387:2206      
    0       RC_pBIN-pROK2   3903    2       39M61S  *       0       0 
    ATGGAGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGGTTTGGGGGGGGTGG    
    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
    AS:i:54 XS:i:50 XN:i:0  XM:i:3  XO:i:0  XG:i:0  NM:i:3  MD:Z:4T4G8G20   YT:Z:UU
    }), 0, 0);
    is_deeply($sam->mismatch_tokens, [
          [ 'M', '4' ],
          [ 'C', 'T' ],
          [ 'M', '4' ],
          [ 'C', 'G' ],
          [ 'M', '8' ],
          [ 'C', 'G' ],
          [ 'M', '20' ],
        ], "forward mismatch_tokens",
    );

    is_deeply($sam->snps, [
            [ 3907, 'T', 'A' ],
            [ 3912, 'G', 'T' ],
            [ 3921, 'G', 'T' ]
        ], "forward snps");
}

#######################################################################
# reverse snps


my %seqlen = (
    CHR1 => 1579921,
    CHR2 => 1579921,
    CHR3 => 1579921,
    CHR4 => 1579921,
    CHR5 => 1579921,
    CHRC => 154478,
    CHRM => 366924,
);

{
    my $sam = Sam::Alignment->new(join("\t", qw{
        chr1:211662:211761:-:211698     
        0       RC_chr1 1368161 255     100M    *       0       0       
AGTTATTAGATAGGTGATTGAGAAAGTGTATATAAATAAATTATTTTAAGAAATTGGTGGAAGCAGGGTTGTAAGGATTTTTTATAGAAATAGTAAGTTG
AGTTATTAGATAGGTGATTGAGAAAGTGTATATAAACAAATCACCTTAAGAAATCGGCGGAAGCAGGGCCGTAAGGACCTTCTACAGAAACAGCAAGTCG    
        XA:i:1  MD:Z:63T36      NM:i:1
        }), \%seqlen, 1);

# Read:   CAACTTACTATTTCTATAAAAAATCCTTACAACCCTGCTTCCACCAATTTCTTAAAATAATTTATTTATATACACTTTCTCAATCACCTATCTAATAACT
# Ref:    CGACTTGCTGTTTCTGTAGAAGGTCCTTACGGCCCTGCTTCCGCCGATTTCTTAAGGTGATTTGTTTATATACACTTTCTCAATCACCTATCTAATAACT
# c2tref:
#                                       ----+----6----+----5----+----4----+----3----+----2----+----1----+----0
#         0----+----1----+----2----+----3----+----4----+----5----+----6

    is($sam->seqid, 'chr1', "reverse seqid");
    is($sam->leftmost, 211662, "reverse leftmost");
    is($sam->rightmost, 211761, "reverse rightmost");

    is_deeply($sam->mismatch_tokens, [
          [ 'M', '36' ],
          [ 'C', 'A' ],
          [ 'M', '63' ],
        ], "reverse mismatch_tokens",
    );

    # say Dumper $sam->snps;

    my ($pos, $base_in_ref, $base_in_read) = @{$sam->snps->[0]};
    # say $sam->readseq;
    is($pos, 211698, "mismatch at correct position");
    is($base_in_ref, 'A');  # rc('T')
    is($base_in_read, 'G'); # rc('C')
}

{
    my $ref = setup_reference(undef, 1);
    my $reads = "$ref.reads";
    my $methsites = "$reads.methsites.gff";
    my $samfile = "$ref.reads.aligned.sam";
    system("perl genome_shear.pl -r .1 -l 100 -n 100 -o $reads $ref");
    system("bowtie -B 1 -S -v 3 $ref.c2t $reads $samfile");

    # load sites from genome_shear.pl
    my %sites;
    my $p = GFF::Parser->new(file => $methsites);
    while (defined(my $gff = $p->next)){
        $sites{uc($gff->sequence), $gff->start}++;
    }

    my $samio = io($samfile)->chomp;
    while (defined(my $line = $samio->getline())){
        next if $line =~ /^@/;
        my $sam = Sam::Alignment->new($line, \%seqlen, 1);
        my $seqid = $sam->seqid;
        for my $snp (@{$sam->snps}) {
            my ($pos, $base_in_ref, $base_in_read) = @$snp;
            ok(exists $sites{uc($seqid), $pos}, "$seqid, $pos methylation where it's supposed to be");
            ok(
                ($base_in_ref eq 'T' && $base_in_read eq 'C')
                ||
                ($base_in_ref eq 'A' && $base_in_read eq 'G'),
                "$seqid, $pos bases indicated methylation (base_in_ref = $base_in_ref, base_in_read = $base_in_read)"
            )
        }
    }
}

done_testing();

