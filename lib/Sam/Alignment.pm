package Sam::Alignment;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Scalar::Util qw/looks_like_number/;

sub new {
    my $class = shift;
    my %opt = @_;

    my $self = bless {}, $class;

    $self->{readid}   = delete $opt{readid};
    $self->{flag}     = delete $opt{flag};
    $self->{seqid}    = delete $opt{seqid};
    $self->{leftmost} = delete $opt{leftmost};
    $self->{mapq}     = delete $opt{mapq};
    $self->{cigar}    = delete $opt{cigar};
    $self->{rnext}    = delete $opt{rnext};
    $self->{pnext}    = delete $opt{pnext};
    $self->{tlen}     = delete $opt{tlen};
    $self->{readseq}  = delete $opt{readseq};
    $self->{readqual} = delete $opt{readqual};
    $self->{optional} = delete $opt{optional};

    return $self;
}

sub new_from_line{
    my $class = shift;
    my $line = shift;

    my ($readid, $flag, $seqid, $leftmost, $mapq, $cigar, $rnext, $pnext, $tlen, $readseq, $readqual, @optional) 
    = split /\t/, $line;

    return $class->new(
        readid   => $readid,
        flag     => $flag,
        seqid    => $seqid,
        leftmost => $leftmost,
        mapq     => $mapq,
        cigar    => $cigar,
        rnext    => $rnext,
        pnext    => $pnext,
        tlen     => $tlen,
        readseq  => $readseq,
        readqual => $readqual,
        optional => \@optional,
    );
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

#######################################################################
# accessors

# 11 mandatory fields are always there
sub readid     { return $_[0]->{readid} }
sub flag       { return $_[0]->{flag} }
sub seqid      { return $_[0]->{seqid} }
sub leftmost   { return $_[0]->{leftmost} }
sub mapq       { return $_[0]->{mapq} }
sub cigar      { return $_[0]->{cigar} }
sub rnext      { return $_[0]->{rnext} }
sub pnext      { return $_[0]->{pnext} }
sub tlen       { return $_[0]->{tlen} }
sub readseq    { return $_[0]->{readseq} }
sub readqual   { return $_[0]->{readqual} }

# two alternatives for leftmost/rightmost when we know 
# the alignment has no gaps, like for bowtie1.
sub dumbstart { return $_[0]->{leftmost} } 
sub dumbend   { return $_[0]->{leftmost} + $_[0]->readlength() } 

# other fields are calculated and maybe memo'd

sub readlength { return length($_[0]->{readseq}) }
sub reverse    { $_[0]->{flag} & $flag_bits{reverse}; }
sub failed_qc  { $_[0]->{flag} & $flag_bits{failed_qc}; }
sub mapped     { !( $_[0]->{flag} & $flag_bits{unmapped} ) }

sub rightmost{
    my $self = shift;

    if (defined $self->{rightmost}){
        return $self->{rightmost};
    }

    my $rightmost = $self->{leftmost};
    my $cigar = $self->{cigar};

    my $length;

    while ($cigar =~ /(\d+)([MIDNSHP=X])/g){
        my $count = $1;
        my $type = $2;
        if ($type ~~ [qw{M D = X N}]){
            $rightmost += $count;
        }
    }
    $rightmost -= 1;
    $self->{rightmost} = $rightmost;

    return $self->{rightmost};
}

# for debugging.  supposed to be equal to readlength
sub _cigarlength{
    my $self = shift;
    my $length;
    my $cigar = $self->{cigar};

    my $cigar_length = 0;

    while ($cigar =~ /(\d+)([MIDNSHP=X])/g){
        my $count = $1;
        my $type = $2;
        if ($type ~~ [qw{M I = X S}]){
            $cigar_length += $count;
        }
    }
    return $cigar_length;

    #croak "BUG: length mismatch between what CIGAR tells us ($cigar_length) and length of SEQ (" . length($seq) . ")" 
    #if length($seq) != $cigar_length;
}

sub mismatches { 
    my $self = shift; 
    if (defined $self->{mismatches}){
        return $self->{mismatches};
    }
    # return $self->{_mismatches} // $self->calc_mismatches() 

    for my $optfield (@{$self->{optional}}) {
        if ($optfield =~ /MD:Z:([A-Z\d\^]+)/){
            # found right optional field.

            my $mdstring = $1;
            my $in_deletion = 0;
            my $position = $self->{leftmost};

            while ($mdstring =~ m{
                    (
                      \d+
                      |
                      \^
                      |
                      [A-Z]
                    )
                }xmg){
                my $token = $1;
                if ($token eq '^'){
                    $in_deletion = 1;
                }
                elsif (looks_like_number $token){
                    $in_deletion = 0;
                    $position += $token;
                }
                # work-in-progress....
                # elsif ( 

                else{
                    $in_deletion = 0;
                }
            }
        }
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
