package BowtieParser;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Moose;
use Carp;
use autodie;    
use open IO  => ":crlf"; # for ARGV...

extends 'Parser';

has base => (
    is => 'ro',
    default => 1,
);

# my ($readid, $strand, $chr, $pos, $read, $quality, $mystery_quantity_that_nobody_understands, $mismatch_string, [PARSED_MATCHES])
# PARSED_MATCHES is [ABS_COORD, BASE_IN_REF, BASE_IN_READ] (everything w.r.t. positive strand).
# all coordinates are base-1

# while (defined(my $alignment = $p->next(1))){
#     my ($readid, $strand, $seqid, $start, $read, $quality, $mystery, $mismatches) = @$alignment;
#     for my $mm (@$mismatches) {
#         my ($abs_coord, $base_in_ref, $base_in_read) = @$mm;
#     }
# }
sub next{
    my $self = shift;
    my $parse_mismatches = shift;
    my $base = $self->base;
    if (defined (my $line = scalar readline $self->filehandle)){
        chomp $line;

        # SOLEXA2_0311_FC62W5LAAXX:2:1:3572:1179#0/1	
        # -	
        # chr5	
        # 365169	
        # TCTTGTNTCNNTATCGGGGGATGCTTGGAGCAGTGACTCC	
        # BBBBBBBBBBBBBB_hhhhhhhefhghhhhehhhgfffff	
        # 0	
        # 29:C>N,30:C>N,33:C>N

        my ($readid, $strand, $seqid, $start, $read, $quality, $mystery, $mismatches_string) = split /\t/, $line, 8;

        $start += (1 - $base); # convert from base-0 to base-1

        if ($parse_mismatches){
            my $rc = $strand eq '+' ? 0 : 1;
            my $len = length($read);

            my @mismatches = map { 
                if (/(\d+):(\w)>(\w)/){
                    # 0 based from 5' end. relative offset can be from 0 to $len - 1
                    my $relative_offset = $1; 
                    my $absolute_offset = $rc ? ($len - $relative_offset - 1): $relative_offset;

                    my $in_ref = $2;
                    my $in_read = $3;

                    # coordinates of mismatches are absolute
                    ([$absolute_offset + $start, $in_ref, $in_read]);
                }
                else{
                    ();
                }
            } split /,/, $mismatches_string;
            return [$readid, $strand, $seqid, $start, $read, $quality, $mystery, \@mismatches];
        }
        else{
            return [$readid, $strand, $seqid, $start, $read, $quality, $mystery];
        }
    }
    return;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;

# SOLEXA2_0311_FC62W5LAAXX:2:1:1584:1185#0/1	-	chr4	14819913	TTGAGGGTTCCAGAGGATTGAGTAAGTGTGGAAATTGGCG	^[`]d^[aKefWbbdZcd]fdWdfacd_]cfffcfccdff	0	
# SOLEXA2_0311_FC62W5LAAXX:2:1:1658:1187#0/1	-	chr5	5937094	TACCTCTTTGCTGGTGTTGTTGATGGAAGGAACATCTGGG	[P^`Ua^^WOfI[eceafae`][`aY\\YSYU^TV\NIU^	0	
# SOLEXA2_0311_FC62W5LAAXX:2:1:2108:1187#0/1	+	chr3	88582	CTGCTTTTTCGCTCGCTGCCACTCTCATCCTCGCCGCGGC	gdggfgggggcgcgfgcgggggagfgcggggg[fgfffab	0	37:T>G,39:T>C
# SOLEXA2_0311_FC62W5LAAXX:2:1:2332:1184#0/1	+	chr1	16710621	GGAGAATGTACAACTTTGTCTCACAGAAGAAACGTAAATA	\]V\LY_aXIa\\^afcfcc`dIa`d^^cO[UP^BBBBBB	0	
# SOLEXA2_0311_FC62W5LAAXX:2:1:2399:1186#0/1	+	chr4	1046365	CTAATCTCCTCTTAAAGACCTTGACAAAATTACATGTAAT	fbccb_ffdfffffdf\affcfaff_dacdeLcRdee[ae	0	
# SOLEXA2_0311_FC62W5LAAXX:2:1:2787:1180#0/1	-	chr1	10472728	AAGATANCTNNAGTGGCCCAAATGGCCAAAATGCTCTGAG	BBBBBBBBBBBB^\agg_fgggggggadgggfggggdggg	0	29:G>N,30:T>N,33:A>N
# SOLEXA2_0311_FC62W5LAAXX:2:1:2884:1185#0/1	-	chr1	7445340	GATTTACTGTAGTAGGTTGTAACAGTTACGCATTTCTCCG	gghghggehgfhhhhhhhghhhhghfhhhhhdghhfffff	0	
# SOLEXA2_0311_FC62W5LAAXX:2:1:3263:1185#0/1	-	chr1	10737588	CCGTGATGAGCAGATTGGCGAGGAACAGAATGCAACACCC	de_ffaa^]Wa[]ffcddWc\[afddcdad_cc[acddad	0	
# SOLEXA2_0311_FC62W5LAAXX:2:1:3378:1182#0/1	+	chr1	10494073	GGAGAATTGTCCCTGGAGGAGAACAGATGNNTANAGCTTG	geggagggfgggeggg_ggagggggBBBBBBBBBBBBBBB	0	29:A>N,30:A>N,33:C>N
# SOLEXA2_0311_FC62W5LAAXX:2:1:3498:1188#0/1	+	chr4	1110659	CTTTTAGAAATTTCTCACTGAACCTTCTCGGAAGCTGGCC	hhhhhhhhhhhhhhhhhhhhgghhghhhhghhgghghhhh	0	
# SOLEXA2_0311_FC62W5LAAXX:2:1:3572:1179#0/1	-	chr5	365169	TCTTGTNTCNNTATCGGGGGATGCTTGGAGCAGTGACTCC	BBBBBBBBBBBBBB_hhhhhhhefhghhhhehhhgfffff	0	29:C>N,30:C>N,33:C>N
# SOLEXA2_0311_FC62W5LAAXX:2:1:3748:1183#0/1	+	chr5	26469425	AGCAGACCCCCTTAACAATCTCTCCGCCATCATTTGTGTA	hhhhhhhhhhhhhhhhhhhhghhfhahheghfhhhfghgh	0	
# SOLEXA2_0311_FC62W5LAAXX:2:1:3902:1181#0/1	-	chr5	3284673	ATCGGANGTNNAAAGGGTTAAATTGTCCCTCAGCAACTGC	[\^^^XB\WBBdeddhhghhhhhhhehhhhhhhhhhhhhh	0	29:C>N,30:C>N,33:G>N
# SOLEXA2_0311_FC62W5LAAXX:2:1:4297:1179#0/1	+	chr2	5546517	GAGAGGAATGTGGTGCTGAAGATAAGTAGNNGCNGCTGCA	hdhgghfhghhghhhhhghhhhgeh`[`BBBBBBBBBBBB	0	29:A>N,30:G>N,33:A>N
# SOLEXA2_0311_FC62W5LAAXX:2:1:4373:1178#0/1	+	chr5	16202278	TGTGGTCAGATAGTGATTTCCTTGCCGGTNNTTNCAATGT	hhghhhhhhhhfhdhhhhhhghhhh^aa[BB_aB\\Z]]Z	0	29:T>N,30:C>N,33:T>N
# SOLEXA2_0311_FC62W5LAAXX:2:1:4439:1188#0/1	-	chr3	19893071	GTGAAATATGAACTCGACAAAACTACGGGTCTCATTAAGG	hhhhhhghhghhh[hh`hghhhhhhhhhhffghhhhhhhh	0	
# SOLEXA2_0311_FC62W5LAAXX:2:1:4581:1184#0/1	+	chr3	2947775	CTTTGATCCTGAAACGAGAGGCCATGTGGAAAAGCTCAAG	fffdf_fcffd^ffageggggfgggf_ffdffafffff]f	0	
# SOLEXA2_0311_FC62W5LAAXX:2:1:4653:1187#0/1	-	chr1	22443079	TCTCCCATGCTAGCATGCACTTGCTTCACCGTCTTATCAT	^aaadd_aa_ac^cKcac^_][dd]a\aaaS[YcYddddY	0	
# SOLEXA2_0311_FC62W5LAAXX:2:1:4866:1186#0/1	-	chr4	14404473	TCTCGGCTGGCTAGAGTGGCATGCATTGCAATGAATCCGC	ddWecggghfdhhghhdghhfhhhghhhhghhhghdfhhf	0	0:G>C,1:T>G
# SOLEXA2_0311_FC62W5LAAXX:2:1:4904:1180#0/1	-	chr2	9364492	AGTTCCNCGNNACAATCTCCTCCGTCACCGTCACCGCAAA	BBBBBBBBBBBBBBBcfcffcfffcfcfffhhghhghhhh	0	29:A>N,30:T>N,33:C>N
# SOLEXA2_0311_FC62W5LAAXX:2:1:5111:1179#0/1	+	chr5	26167329	GAGAAGATTCCTCTCTTCTACGGTCAGCTNNCCNGCCGGT	hahhdhhhhhfhghhhhhhhhhghh]Z]ZBB_XBXXW\XM	0	29:C>N,30:A>N,33:G>N
# SOLEXA2_0311_FC62W5LAAXX:2:1:5248:1179#0/1	-	chr3	2410777	CTTAGGNTGNNCAGAACGGTTGGTGATGAATCCATTGAAG	BBBBBBBBBBBBBBBcagggffffffdfffaggggggggg	0	29:G>N,30:T>N,33:A>N
# SOLEXA2_0311_FC62W5LAAXX:2:1:5281:1186#0/1	-	chr1	24220866	TCCACAACCAACCCATCTCTAAACCTAAACTGTTTGTTTT	_ac^caffffbfefcdebdd_dfffdffff]fffdf_cff	0	4:T>G
# SOLEXA2_0311_FC62W5LAAXX:2:1:5440:1180#0/1	+	chr3	16537602	CCCCTCTTCAGAGCAACGCTACGGGGGAGNNTTNTATGTA	\fdgcfcfffg[gfgcfdffcfcacBBBBBBBBBBBBBBB	0	29:A>N,30:G>N,33:C>N
# SOLEXA2_0311_FC62W5LAAXX:2:1:5486:1188#0/1	-	chr5	2792548	TTCTGGATTCGCAACAAGATTGATGGGTGTGAAAGTCCTT	_ccacegfggafgggggcggegggffffffcgggggggef	0	
# SOLEXA2_0311_FC62W5LAAXX:2:1:5624:1188#0/1	+	chr3	20945227	CGAGAATCTGAAGGTAATTTGCGTAATTCAAGCTCGTAGC	_`Z`WVYUV\OQ]VM_V_\_fWff_[^\YY]`YT`BBBBB	0	
# SOLEXA2_0311_FC62W5LAAXX:2:1:5761:1185#0/1	+	chr3	6476870	GTTTTATGAATGGGATCAAAGTTTCTTTTTTTCTTTTATA	hhhhhfcfcdffffdfheehghfgghghhhhhbhhfc`c]	0	
# SOLEXA2_0311_FC62W5LAAXX:2:1:5878:1186#0/1	-	chr5	26413865	CGTTCAGTCTATCGGTATCTAACTCCATGCCAAGGCCATC	hehhhhgefghhhfhhfhhhhhhghhhhhghhhhhghhhh	0	
# SOLEXA2_0311_FC62W5LAAXX:2:1:5920:1182#0/1	-	chr2	2359495	CATTAGCTTAGCTTCTTCCTCATCACTCCCTTGCCTTCCC	cggggefcfffd[fffcff[ee`e`eWceefafcff[ffd	0	0:A>C
# SOLEXA2_0311_FC62W5LAAXX:2:1:6028:1188#0/1	+	chr2	14529683	GGAGGTTCGATTTCTATACAAAGATACCACAAAGGAAAGT	hhfhhhhhghhhhhhhfhhhhhhhhhhhghgdhhfhhhgf	0	0:A>G
