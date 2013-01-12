package Sam::Alignment;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Scalar::Util qw/looks_like_number/;
use Mouse;

around BUILDARGS => sub{
    my ($orig, $class, $line) = @_;
    chomp($line);

    my ($readid, $flag, $seqid, $leftmost, $mapq, $cigar, $rnext, $pnext, $tlen, $readseq, $readqual, @optional) 
    = split /\t/, $line;

    return $class->$orig(
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
};

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
has readid   => ( is => 'ro', required => 1 );
has flag     => ( is => 'ro', required => 1 );
has seqid    => ( is => 'ro', required => 1 );
has leftmost => ( is => 'ro', required => 1 );
has mapq     => ( is => 'ro', required => 1 );
has cigar    => ( is => 'ro', required => 1 );
has rnext    => ( is => 'ro', required => 1 );
has pnext    => ( is => 'ro', required => 1 );
has tlen     => ( is => 'ro', required => 1 );
has readseq  => ( is => 'ro', required => 1 );
has readqual => ( is => 'ro', required => 1 );
has optional => ( is => 'ro', required => 1 );

# two alternatives for leftmost/rightmost when we know 
# the alignment has no gaps, like for bowtie1.
sub dumbstart { return $_[0]->leftmost } 
sub dumbend   { return $_[0]->leftmost + $_[0]->readlength() } 

# other fields are calculated and maybe memo'd

sub readlength { return length($_[0]->readseq) }
sub reverse    { $_[0]->flag & $flag_bits{reverse}; }
sub failed_qc  { $_[0]->flag & $flag_bits{failed_qc}; }
sub mapped     { !( $_[0]->flag & $flag_bits{unmapped} ) }

has rightmost => ( is => 'ro', lazy_build => 1);

sub _build_rightmost{
    my $self = shift;

    my $rightmost = $self->leftmost;
    my $cigar = $self->cigar;

    my $length;

    while ($cigar =~ /(\d+)([MIDNSHP=X])/g){
        my $count = $1;
        my $type = $2;
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
    my $cigar = $self->cigar;
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

has mismatch_string => ( is => 'ro', lazy_build => 1,);
 
sub _build_mismatch_string{
    my $self = shift;
    for my $optfield (@{$self->optional}) {
        if ($optfield =~ /MD:Z:([A-Z\d\^]+)/){
            $self->{mismatch_string} = $1;
            return $1;
        }
    }
    return;
}

# sub mismatches { 
#     my $self = shift; 
#     if (defined $self->{mismatches}){
#         return $self->{mismatches};
#     }
#     my $mdstring = $self->_mismatch_string;
# 
#     my $mdstring = $1;
#     my $in_deletion = 0;
#     my $position = $self->{leftmost};
# 
#     while ($mdstring =~ m{(
#             \d+
#             |
#             \^
#             |
#             [A-Z]
#             )}xmg){
#         my $token = $1;
#         if ($token eq '^'){
#             $in_deletion = 1;
#         }
#         elsif (looks_like_number $token){
#             $in_deletion = 0;
#             $position += $token;
#         }
#         # work-in-progress....
#         # elsif ( 
# 
#         else{
#             $in_deletion = 0;
#         }
#     }
# }

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
