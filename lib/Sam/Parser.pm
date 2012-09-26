package Sam::Parser;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Hash::Util qw/lock_keys unlock_keys/;
use parent 'ParserNoMoose';
use List::MoreUtils qw/all/;

sub new {
    my $class = shift;
    my %opt = @_;

    my $convert_rc = delete $opt{convert_rc};
    my $skip_unmapped = delete $opt{skip_unmapped};

    my $self = $class->SUPER::new(%opt);

    unlock_keys(%$self);

    $self->{length}          = {};
    $self->{rc_sequence}     = {};    # rc_sequence caches the names of the ^RC_ chromosomes
    $self->{program_name}    = undef;
    $self->{program_version} = undef;
    $self->{command_line}    = undef;
    $self->{sort_order}      = undef;
    $self->{sam_version}     = undef;
    $self->{convert_rc}      = $convert_rc // 0;
    $self->{skip_unmapped}   = $skip_unmapped // 1;
    $self->{putback}         = undef; # during constructor, read until first alignment and putback here
                                      # want to use seek(), but doesn't work on stdin properly (?)

    lock_keys(%$self);

    # read headers, putback first alignment line into $self->{putback}
    HEADER:
    while (defined(my $line = readline $self->{handle})){
        chomp $line;
        if ($line =~ /^@/){
            $self->parse_header($line);
        }
        else{
            $self->{putback} = $line;
            last HEADER;
        }
    }

    return $self;
}

sub next{
    my $self = shift;

    # process putback
    if (defined (my $pb = $self->{putback})){
        undef $self->{putback};
        my $align = $self->parse_alignment($pb);
        if ($align->{mapped} || ! $self->{skip_unmapped}){
            return $align;
        }
    }

    while (defined(my $line = readline $self->{handle})){
        chomp $line;
        my $align = $self->parse_alignment($line);
        if ($align->{mapped} || ! $self->{skip_unmapped}){
            return $align;
        }
    }

    return;
}

# @HD	VN:1.0	SO:unsorted
# @SQ	SN:chr1	LN:30432563
# @SQ	SN:chr2	LN:19705359
# @SQ	SN:chr3	LN:23470805
# @SQ	SN:chr4	LN:18585042
# @SQ	SN:chr5	LN:26992728
# @SQ	SN:chrc	LN:154478
# @SQ	SN:chrm	LN:366924
# @PG	ID:Bowtie	VN:0.12.7	CL:"bowtie -S -f -B 1 -v 3 --best /wip/tools/genomes/AT/TAIR_reference.fas read.rc.fasta"
sub parse_header{
    my $self = shift;
    my $line = shift;

    my ($type_string, @parts) = split /\s+/, $line;

    my ($type) = $type_string =~ /^@(\w\w)/;

    die "can't parse header type" unless $type;

    my %header;

    for my $part (@parts) {
        if ($part =~ /(\w\w):(.*)/){
            my ($key, $value) = ($1, $2);
            $header{$key} = $value;
        }
    }

    if ($type eq 'HD' and defined($header{VN})){
        $self->{sam_version} = $header{VN};
        if ($header{SO}){
            $self->{sort_order} = $header{SO};
        }
    }
    elsif ($type eq 'SQ' and all { defined $header{ $_ } } qw/SN LN/){
        $self->{length}{$header{SN}} = $header{LN};
        if ($self->{convert_rc} and $header{SN} =~ /^RC_/){
            $self->{rc_sequence}{$header{SN}} = 1;
        }
    }
    elsif ($type eq 'RG' and defined($header{ID})){
        # not sure what this is but leaving as stub
    }
    elsif ($type eq 'PG' and defined($header{ID})){
        $self->{program_name} = $header{ID};
        if ($header{CL}){
            $self->{command_line} = $header{CL};
        }
        if ($header{VN}){
            $self->{program_version} = $header{VN};
        }
    }
    elsif ($type eq 'CO'){
        # ignore
    }
    else{
        die "can't parse header line";
    }
}

our %flag_bits = (
    multiple_segments   => 0x1,
    each_segment_aligns => 0x2,
    unmapped            => 0x4,
    next_unmapped       => 0x8,
    reverse             => 0x10,
    next_reverse        => 0x20,
    first_segment       => 0x40,
    last_segment        => 0x80,
    secondary_alignment => 0x100,
    failed_qc           => 0x200,
    duplicate           => 0x400,
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

# return line as-is for header lines. otherwise, hashref
sub parse_alignment{
    my $self = shift;
    my $line = shift;
    $line =~ tr/\n\r//d;

    my ($qname, $flag_number, $rname, $leftmost, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual) 
    = split /\t/, $line;

    my $rightmost;
    my $length;

    # parse flag
    my $reverse   = $flag_number & $flag_bits{reverse};
    my $mapped    = !($flag_number & $flag_bits{unmapped});
    my $failed_qc = $flag_number & $flag_bits{failed_qc};


    if (! $mapped){
        undef($leftmost);
        undef($rname);
    }
    else{

        # determine the rightmost and length (in genome)
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


    my $samline = bless({
            mapped    => $mapped,

            seqid     => $rname,
            length    => $length,
            leftmost  => $leftmost,
            rightmost => $rightmost,

            reverse   => $reverse,
            readid    => $qname,
            seq       => $seq,
            cigar     => $cigar,
            mapq      => $mapq,
            failed_qc => $failed_qc,

        }, "SAMLINE");

    if (exists $self->{rc_sequence}{$qname}){
        return $self->reverse_sam($samline);
    }
    else{
        return $samline;
    }
}

#######################################################################
# reversal code

sub reverse_sam{
    # rname: remove RC_
    # leftmost/rightmost: swap and convert
    # cigar: reverse_cigar
    # seq,qual: reverse complement
    # pnext: not sure
    # rnext: not sure
    # opt/MD: like cigar
}

sub reverse_cigar{
    #my ($self, $seq, $cig) = @_;
}

sub reverse_md{
    #my ($self, $seq, $cig) = @_;
}

1;

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

