package Sam::Alignment;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Scalar::Util qw/looks_like_number/;
use Mouse;

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

around BUILDARGS => sub{
    my ($orig, $class, $line, $seqlengths, $tryfixrc) = @_;
    chomp($line);

    my ($readid, $flag, $seqid, $leftmost, $mapq, $cigar, $rnext, $pnext, $tlen, $readseq, $readqual, @optional) 
    = split /\t/, $line;

    my $mismatch_string;
    my $edit_distance;
    for my $o (@optional) {
        if ($o =~ /MD:Z:([A-Z\d\^]+)/){
            $mismatch_string = $1;
        }
        elsif ($o =~ /NM:i:(\d+)/){
            $edit_distance = $1;
        }
    }

    croak "If \$tryfixrc is true, must also pass \$seqlengths" if $tryfixrc && ! $seqlengths;

    my $fixrc;
    if ($tryfixrc && $seqid =~ s/^RC_//){
        $fixrc = 1;
        croak "\$tryfixrc is true, but \$seqid $seqid not found in \$seqlengths?"
        if (!exists $seqlengths->{$seqid});

        # $readid unchanged
        $flag = $flag & ~$flag_bits{reverse};
        # $seqid already stripped of RC_
        $leftmost = $seqlengths->{$seqid} - $leftmost + 1;
        # $mapq unchanged
        # $cigar should be reversed, but should be done lazily since we probably won't need it
        # not sure about rnext 
        # not sure about pnext 
        # tlen probably unchanged 
        $readseq = reverse_complement($readseq);
        $readqual = scalar reverse $readseq;
        
        # edit_distance can remain same
        # $mismatch_string should be reversed like cigar lazily.
    }

    return $class->$orig(
        readid   => $readid,
        flag     => $flag,
        seqid    => $seqid,
        leftmost => $leftmost,
        mapq     => $mapq,
        original_cigar_string => $cigar,
        rnext    => $rnext,
        pnext    => $pnext,
        tlen     => $tlen,
        readseq  => $readseq,
        readqual => $readqual,

        # from optional
        original_mismatch_string => $mismatch_string,
        edit_distance            => $edit_distance,

        fixrc    => $fixrc,
    );
};


#######################################################################
# accessors

# 11 mandatory fields are always there
has readid   => ( is => 'ro', required => 1 );
has flag     => ( is => 'ro', required => 1 );
has seqid    => ( is => 'ro', required => 1 );
has leftmost => ( is => 'ro', required => 1 );
has mapq     => ( is => 'ro', required => 1 );
has original_cigar_string    => ( is => 'ro', required => 1 );
has rnext    => ( is => 'ro', required => 1 );
has pnext    => ( is => 'ro', required => 1 );
has tlen     => ( is => 'ro', required => 1 );
has readseq  => ( is => 'ro', required => 1 );
has readqual => ( is => 'ro', required => 1 );

has original_mismatch_string => ( is => 'ro', required => 1,);
has edit_distance   => ( is => 'ro', required => 1,);

# true if tryfixrc and /^RC_/
has fixrc => ( is => 'ro', required => 1 );

# two alternatives for leftmost/rightmost when we know 
# the alignment has no gaps, like for bowtie1.
sub dumbstart { return $_[0]->leftmost } 
sub dumbend   { return $_[0]->leftmost + $_[0]->readlength() } 

# other fields are calculated and maybe memo'd

sub readlength { return length($_[0]->readseq) }
sub reverse    {    $_[0]->flag & $flag_bits{reverse};     }
sub failed_qc  {    $_[0]->flag & $flag_bits{failed_qc};   }
sub mapped     { !( $_[0]->flag & $flag_bits{unmapped}   ) }

has rightmost => ( is => 'ro', lazy_build => 1);

sub _build_rightmost{
    my $self = shift;

    my $rightmost = $self->leftmost;
    my $length;

    for my $c (@{$self->cigar}) {
        my ($count, $type) = @$c;
        if ($type ~~ [qw{M D = X N}]){
            $rightmost += $count;
        }
    }
    $rightmost -= 1;

    return $rightmost;
}

# for debugging.  supposed to be equal to read length
has cigarlength => ( is => 'ro', lazy_build => 1,);

sub _build_cigarlength{
    my $self = shift;
    my $cigar_length = 0;
    for my $c (@{$self->cigar}) {
        my ($count, $type) = @$c;

        if ($type ~~ [qw{M I = X S}]){
            $cigar_length += $count;
        }
    }
    return $cigar_length;

    #croak "BUG: length mismatch between what CIGAR tells us ($cigar_length) and length of SEQ (" . length($seq) . ")" 
    #if length($seq) != $cigar_length;
}

has cigar => ( is => 'ro', lazy_build => 1 );

sub _build_cigar{
    my $self = shift;
    my $cigar = $self->original_cigar_string;

    return if $cigar eq '*';

    my @accum;

    while ($cigar =~ /(\d+)([MIDNSHP=X])/g){
        my $count = $1;
        my $type = $2;
        push @accum, [$count, $type];
    }
    if ($self->fixrc){
        # reverse accum?
    }
    return \@accum;
}

has mismatches => ( is => 'ro', lazy_build => 1 );

sub _build_mismatches { 
    my $self = shift; 

    my $mdstring = $self->original_mismatch_string;

    my $in_deletion = 0;
    my $position = $self->{leftmost};

    my @accum;

    # M = match
    # C = change
    # D = change

    while ($mdstring =~ m{( \d+ | \^ | [A-Z] )}xmg){
        my $token = $1;
        if ($token eq '^'){
            $in_deletion = 1;
        }
        elsif (looks_like_number $token){
            $in_deletion = 0;
            $position += $token;

            push @accum, ['M', $token];
        }
        else{ # [A-Z]
            if ($in_deletion){
                push @accum, ['D', $token];
                # no in_deletion reset, keep going
            }
            else{
                push @accum, ['C', $token];
            }
        }
    }

    if ($self->fixrc){
        # reverse accum?
    }
    return \@accum;
}

#######################################################################

# sub reverse_cigar{
#     #my ($self, $seq, $cig) = @_;
# }
# 
# sub reverse_md{
#     #my ($self, $seq, $cig) = @_;
# }


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

#######################################################################
# stringification


use overload '""' => \&stringify;

sub stringify{
    my $self = shift;
    my $mapped    = $self->mapped;
    my $rname     = $self->seqid;
    my $length    = $self->length;
    my $leftmost  = $self->leftmost;
    my $rightmost = $self->rightmost;
    my $failed_qc = $self->failed_qc;
    my $reverse   = $self->reverse;
    my $qname     = $self->readid;
    my $seq       = $self->seq;
    my $cigar     = $self->cigar;
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

=head1 Sam::Alignment

=head2 $sam->readid 

Read ID from FASTQ. (Column 1)

=head2 $sam->flag 

Bitfield of flags. (Column 2)

=head2 $sam->mapped

Did the read align? Extracted from flags.

=head2 $sam->reverse

Was the alignment on the reverse strand? Note that readseq and positions are
always w.r.t. to 5' end. Extracted from flags.

=head2 $sam->failed_qc

Does the read have high enough quality? Extracted from flags.

=head2 $sam->seqid 

Chromosome name from reference. (Column 3)

=head2 $sam->leftmost 

Leftmost position. (Column 4)

=head2 $sam->rightmost 

Rightmost position.  Calculated from leftmost and CIGAR string.

=head2 $sam->length 

Length. Rightmost-leftmost+1.

=head2 $sam->mapq 

Mapping quality. (Column 5)

=head2 $sam->cigar 

See CIGAR section. (Column 6)

=head2 $sam->rnext 

"Ref. name of the mate/next segment". Unused. (Column 7)

=head2 $sam->pnext 

"Position of the mate/next segment". Unused. (Column 8)

=head2 $sam->tlen 

"observed Template LENgth". Unused. (Column 9)

=head2 $sam->readseq 

Actual read sequence from FASTQ. (ACGTN's) (Column 10)

=head2 $sam->readqual 

Actual read quality from FASTQ. (Solexa/etc codes) (Column 11)

=head2 $sam->mismatches

Mismatch position.  Extracted from MD:Z: field of column 12.


=cut
