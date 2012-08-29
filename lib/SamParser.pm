package SamParser;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(parse_sam_line);

# not used, since only need a few...
my %flags = (
    0x1 => "multiple_segments",
    0x2 => "each_segment_aligns",
    0x4 => "unmapped",
    0x8 => "next_unmapped",
    0x10 => "reverse",
    0x20 => "next_reverse",
    0x40 => "first_segment",
    0x80 => "last_segment",
    0x100 => "secondary_alignment",
    0x200 => "failed_qc",
    0x400 => "duplicate",
);

# sample output:
# $VAR1 = {
#           'readid' => 'SALK_027726.22.20.x',
#           'cigar' => '93M1D107M1I66M1I54M',
#           'failed_qc' => 0,
#           'seqid' => 'chr5',
#           'reverse' => 16,
#           'unmapped' => 0,
#           'length' => 322,
#           'leftmost' => '25847254',
#           'seq' => 'AGCTTAGTCCATTCATCCAATCTTCTATAAACTGCAGTTTCACAAGAGTTGTGATACTCTATCCAACTTGAATTGCAAAGCCCTCTAATTCGATCTTCGAAATTTCAGCTCCGATCTATAGGCCCTAGATTCACAGATTTCAGATAGAGAAAGGAGGGAACGACATGCTCTCTGTGTAATGACCGTAGTACTCGTAATCAGTTGTGACCACCATAATAATTCTTTGGGCGCTCGATAACCTTCACCGGATGGTATTTTCCGTATAACACATTGACACATCACGCGAGTGCCAACCTCGAGCGGAAGCACCGCCATTCTTCTC',
#           'rightmost' => 25847574
#         };
#
# return line as-is for header lines. otherwise, hashref
sub parse_sam_line{
    my $line = shift;
    $line =~ tr/\n\r//d;

    if ($line =~ /^@/){
        return $line;
    }
    else{
        my %result;
        
        my ($qname, $flag, $rname, $leftmost, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual) 
        = split /\t/, $line;

        my $unmapped = ($flag & 0x4);
        my $reverse = ($flag & 0x10);
        my $failed_qc = ($flag & 0x200);

        my $rightmost;
        my $length;

        if ($unmapped){
            undef($leftmost);
            undef($rightmost);
            undef($rname);
            undef($length);
        }
        else{
            $rightmost = $leftmost;
            my $cigar_length = 0;

            while ($cigar =~ /(\d+)([MIDNSHP=X])/g){
                my $count = $1;
                my $type = $2;
                if ($type ~~ [qw{M D = X N}]){
                    $rightmost += $count;
                }
                if ($type ~~ [qw{M I = X S}]){
                    $cigar_length += $count;
                }
            }
            $rightmost -= 1;
            croak "BUG: length mismatch between what CIGAR tells us ($cigar_length) and length of SEQ (" . length($seq) . ")" 
            if length($seq) != $cigar_length;

            $length = $cigar_length;
        }

        return bless({
            mapped    => ! $unmapped,

            seqid     => $rname,
            length    => $length,
            leftmost  => $leftmost,
            rightmost => $rightmost,

            failed_qc => $failed_qc,
            reverse   => $reverse,
            readid    => $qname,
            seq       => $seq,
            cigar     => $cigar,
            mapq      => $mapq,
        }, "SAMLINE");
    }
}

package SAMLINE;
use strict;
use warnings FATAL => "all";
use 5.010_000;

use overload '""' => \&stringify;
sub stringify{
    my $self = shift;
    my $mapped    = $self->{mapped};
    my $rname     = $self->{seqid};
    my $length    = $self->{length};
    my $leftmost  = $self->{leftmost};
    my $rightmost = $self->{rightmost};
    my $failed_qc = $self->{failed_qc};
    my $reverse   = $self->{reverse};
    my $qname     = $self->{readid};
    my $seq       = $self->{seq};
    my $cigar     = $self->{cigar};
    my $strand    = $reverse ? "-" : "+";

    my $seqid = $qname eq '0' ? '*' : $qname;

    if ($mapped){
        return "$seqid mapped to $rname from $leftmost to $rightmost on $strand strand ($cigar)";
    }
    else{
        return "$seqid did not map";
    }
}


1;


# CIGAR: Description
# M alignment match (can be a sequence match or mismatch)
# I insertion to the reference
# D deletion from the reference
# N skipped region from the reference
# S soft clipping (clipped sequences present in SEQ)
# H hard clipping (clipped sequences NOT present in SEQ)
# P padding (silent deletion from padded reference)
# = sequence match
# X sequence mismatch

# manual = http://samtools.sourceforge.net/SAM1.pdf

# notes:
# rightmost = POS + #(M/D/=/X/N) - 1
# definitely not not s,h,p... 
# not I b/c that only advances read
# manual says N is intron (for mrna). 
# I think =/X are subtype of M?
#
# length of SEQ = M/I/S/=/X (from manual)

