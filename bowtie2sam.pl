#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use IO::All;
use Pod::Usage;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/lib";
use BowtieParser;
use FastaReader;
use DZUtil qw/reverse_complement/;

my $result = GetOptions(
    "help"          => \my $help,
    "reference|r=s" => \(my $reference),
    "base|b=i"      => \(my $base = 1),
);
if (! $result || !$reference || (0 == @ARGV && -t STDIN)){
    pod2usage(-verbose => 2, -noperldoc => 1);
}

my $fr = FastaReader->new(file => $reference, slurp => 0);
my $bp = BowtieParser->new(file => \*ARGV);

say '@HD	VN:1.0	SO:unsorted';

# @SQ	SN:Chr1	LN:30432563
# @SQ	SN:Chr2	LN:19705359
# @SQ	SN:Chr3	LN:23470805
# @SQ	SN:Chr4	LN:18585042
# @SQ	SN:Chr5	LN:26992728
# @SQ	SN:ChrC	LN:154478
# @SQ	SN:ChrM	LN:366924
my %seqlen = $fr->sequence_lengths();
for my $seq (sort keys %seqlen) {
    say "\@SQ\tSN:$seq\tLN:$seqlen{$seq}";
}

while (defined(my $alignment = $bp->next_fixrc(1, $fr))){
    my ($readid, $strand, $seqid, $start, $read, $quality, $mystery, $mismatches) = @$alignment;
    my $len = length($read);
    my $num_mm = scalar(@$mismatches);

    say join("\t", 
        $readid,          # HS2:306:C1A3MACXX:1:1101:1331:2213#/1
        ($strand eq '-' ? 16 : 0), # flags
        $seqid,           # seqid
        $start,           # 19230337
        255,              # 255 (?)
        "${len}M",        # 50M (cigar)
        "*",              # *
        0,                # 0
        0,                # 0
        $read,            # TATTTTNGTTATTGTATGTGAATTGTTTATTTTGAATTATATGATTTTTT
        $quality,         # IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
        build_mdz_string($mismatches, $len, $strand eq '-'), # MD:Z:6T43, mismatch positions
        "NM:i:$num_mm",  
    );
}
use Scalar::Util qw/looks_like_number/;

sub build_mdz_string{
    my ($mismatches, $len) = @_;
    my $position = 0;
    my @accum;
    for my $mm (@$mismatches) {
        my (undef, $base_in_ref, undef, $rel_coord) = @$mm;
        push @accum, $rel_coord - $position, $base_in_ref;
        $position = $rel_coord + 1;
    }
    if ($position != $len){
        push @accum, $len - $position;
    }
    for my $x (@accum) {
        if (looks_like_number($x) and $x < 0){
            die "$x < 0";
        }
    }
    return "MD:Z:" . join '', @accum;
}

# HS2:306:C1A3MACXX:1:1101:6987:2227#/1   
# +       
# Chr3    
# 13594770        
# AAGAATTTAAATTATAATTTGATTTTGTAAGTTTAAGAAGTGTATTTTTG      
# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      
# 2       
# 2:T>G,37:T>A

# HS2:306:C1A3MACXX:1:1101:6987:2227#/1   
# 0       
# Chr3    
# 13594770        
# 255     
# 50M     
# *       
# 0       
# 0       
# AAGAATTTAAATTATAATTTGATTTTGTAAGTTTAAGAAGTGTATTTTTG      
# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      
# XA:i:2  
# MD:Z:2T34T12    
# NM:i:2

# ==> foo.txt <==
# @HD	VN:1.0	SO:unsorted
# @SQ	SN:Chr1	LN:30432563
# @SQ	SN:Chr2	LN:19705359
# @SQ	SN:Chr3	LN:23470805
# @SQ	SN:Chr4	LN:18585042
# @SQ	SN:Chr5	LN:26992728
# @SQ	SN:ChrC	LN:154478
# @SQ	SN:ChrM	LN:366924
# @SQ	SN:RC_Chr1	LN:30432563
# @SQ	SN:RC_Chr2	LN:19705359
# @SQ	SN:RC_Chr3	LN:23470805
# @SQ	SN:RC_Chr4	LN:18585042
# @SQ	SN:RC_Chr5	LN:26992728
# @SQ	SN:RC_ChrC	LN:154478
# @SQ	SN:RC_ChrM	LN:366924
# @PG	ID:Bowtie	VN:0.12.7	CL:"bowtie /home/toshiro/genomes/AT/TAIR_reference.fas.c2t -f -B 1 -v 2 -5 0 -3 50 --best --strata -k 10 -m 10 --norc -S - foo.txt"
# HS2:306:C1A3MACXX:1:1101:1381:2127#/1	0	ChrC	114937	255	50M	*	0	0	NTTAATAGTTTTTTTAGTGATTATATTTTGTAAAGTTGGATTTATTTTTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	XA:i:1	MD:Z:0A49	NM:i:1
# HS2:306:C1A3MACXX:1:1101:1305:2195#/1	4	*	0	0	*	*	0	0	ATTGANNNNNNGAANTTAAATTGTAATTTGATTTTAAAGGTGTAAGAATT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	XM:i:0
# HS2:306:C1A3MACXX:1:1101:1331:2213#/1	0	RC_Chr3	19230337	255	50M	*	0	0	TATTTTNGTTATTGTATGTGAATTGTTTATTTTGAATTATATGATTTTTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	XA:i:1	MD:Z:6T43	NM:i:1
# HS2:306:C1A3MACXX:1:1101:1407:2238#/1	0	ChrC	143669	255	50M	*	0	0	GATTATTGTGATTGTTTAGGATATTTTTTTTAGTTTTTAGAATTTATTTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	XA:i:0	MD:Z:50	NM:i:0
# 
# ==> test/reads.fastq_1-50.bowtie <==
# HS2:306:C1A3MACXX:1:1101:1381:2127#/1	+	ChrC	114937	NTTAATAGTTTTTTTAGTGATTATATTTTGTAAAGTTGGATTTATTTTTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	0	0:A>N
# HS2:306:C1A3MACXX:1:1101:1331:2213#/1	+	RC_Chr3	19230337	TATTTTNGTTATTGTATGTGAATTGTTTATTTTGAATTATATGATTTTTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	0	6:T>N
# HS2:306:C1A3MACXX:1:1101:1407:2238#/1	+	ChrC	143669	GATTATTGTGATTGTTTAGGATATTTTTTTTAGTTTTTAGAATTTATTTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	1	
# HS2:306:C1A3MACXX:1:1101:1407:2238#/1	+	RC_ChrC	59499	GATTATTGTGATTGTTTAGGATATTTTTTTTAGTTTTTAGAATTTATTTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	1	
# HS2:306:C1A3MACXX:1:1101:1571:2189#/1	+	ChrC	125794	AAGAGTATAGGTTATAATTTATTTTTTTTTTATTGGATATTGAATTTTTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	0	41:A>G
# HS2:306:C1A3MACXX:1:1101:1713:2195#/1	+	RC_Chr5	11348274	TTATATGTTAATTTTTAAGATGATGTAAATATGAAATTAATTTTTATAAG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	0	
# HS2:306:C1A3MACXX:1:1101:1550:2238#/1	+	Chr2	9927	TAGTTTTTAAGTTAAGAAAGTTGTAAAAGTTAAGAATTAGTATTAAATGA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	0	
# HS2:306:C1A3MACXX:1:1101:1828:2087#/1	+	ChrC	89037	NTGAGGTTATTGATAAAAAAGATTTGTTTAAGTTATTTTGTTTTTTTTTG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	1	0:T>N
# HS2:306:C1A3MACXX:1:1101:1828:2087#/1	+	RC_ChrC	4867	NTGAGGTTATTGATAAAAAAGATTTGTTTAAGTTATTTTGTTTTTTTTTG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	1	0:T>N
# HS2:306:C1A3MACXX:1:1101:1908:2109#/1	+	Chr4	339266	NAATATTGTTTTTAGAGTTGTTAGGAGAATTTATTTTTATTAGAGAGTAG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	0	0:A>N
# HS2:306:C1A3MACXX:1:1101:1859:2135#/1	+	Chr3	5339942	NTAAGGTAAGTGTAGGGTTTAAAGTTTTATAATTTTATATGATTATTATA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	0	0:A>N
# HS2:306:C1A3MACXX:1:1101:1845:2173#/1	+	RC_Chr2	19700889	TTAATTTTTTTAAAGTAATAGTGTTGGAGGTATGATTTGGTTAATTAAGA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	1	
# HS2:306:C1A3MACXX:1:1101:1845:2173#/1	+	RC_Chr3	9261378	TTAATTTTTTTAAAGTAATAGTGTTGGAGGTATGATTTGGTTAATTAAGA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	1	
# HS2:306:C1A3MACXX:1:1101:1773:2194#/1	+	Chr5	10264225	GTTGAAAAATTGTTATATATTTTGATTATTAAAATGTTATTTTTGTTATA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	0	
# HS2:306:C1A3MACXX:1:1101:1865:2195#/1	+	ChrC	95518	TTTATGTTTGTTTGAGTAATAGTAATGAGATTTTTGAATATTATGTTAAG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	1	
# HS2:306:C1A3MACXX:1:1101:1865:2195#/1	+	RC_ChrC	11348	TTTATGTTTGTTTGAGTAATAGTAATGAGATTTTTGAATATTATGTTAAG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	1	
# HS2:306:C1A3MACXX:1:1101:2100:2165#/1	+	Chr5	4134059	ATGGTTTTAGATAGATATAGTAAGTATTTGTTATAATTGATATTGTTGTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	0	
# HS2:306:C1A3MACXX:1:1101:2141:2171#/1	+	Chr2	17606203	TGTTGTTGGGGTTGTTAATAATTTATATGGTTAAGTTTAATTGTTATTGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	0	
# HS2:306:C1A3MACXX:1:1101:2022:2228#/1	+	RC_Chr1	12975864	TTAATAATTATTGAGAGATTTGTTAATTTTAAGTATATTAAAAAATAAAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	0	
# HS2:306:C1A3MACXX:1:1101:2320:2113#/1	+	Chr4	8537583	NAATATTTTTTTGTTATGAGGTTTGAGGGTTGATTATGATTTTGTTTTTA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	0	0:A>N

=head1 bowtie2sam.pl 

Usage examples:

 bowtie2sam.pl -r reference.fas input.bowtie > output.sam

=cut
