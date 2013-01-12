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

sub check_basics{
    my $name = shift;
    my $is_mapped = shift;
    my $line = join "\t", @_;

    my $sam = Sam::Alignment->new($line);
    my $p = "$name (basic)"; # prefix
    if ($is_mapped){
        ok($sam->mapped, "$p: should've mapped");
        ok(looks_like_number($sam->leftmost) && $sam->leftmost > 0, "$p: positive position");
        isnt($sam->cigar, '*', "$p: has a cigar string");
        like($sam->cigar, qr/^(?:\d+[MIDNSHP=X])+$/, "$p: cigar looks valid");
        isnt($sam->seqid, '*', "$p: has a seqid");

        is(
            $sam->cigarlength(), 
            $sam->readlength(),
            "$p: cigar length is readlength"
        );

        ok($sam->mismatch_string, "$p: has mismatch string");
    }
    else{
        ok(! $sam->mapped, "$p: shouldn't have mapped");
        is($sam->cigar, '*', "$p: doesn't have a cigar string");
        is($sam->seqid, '*', "$p: doesn't have a seqid");
    }
    is(length($sam->readseq), length($sam->readqual), "$name (basic): readseq and readqual same length");

    return $sam;
}

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


done_testing();
