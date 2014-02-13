package Sam::Alignment;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Scalar::Util qw/looks_like_number/;
use DZUtil qw/reverse_complement/;
use Mouse;

our %flag_bits = (
    multiple_segments   => 0x1,   # for paired
    each_segment_aligns => 0x2,
    unmapped            => 0x4,
    next_unmapped       => 0x8,
    reverse             => 0x10,
    next_reverse        => 0x20,
    first_segment       => 0x40,  # for paired
    last_segment        => 0x80,  # for paired
    secondary_alignment => 0x100,
    failed_qc           => 0x200,
    duplicate           => 0x400,
);

# my $sam = Sam::Alignment->new($line, \%sequence_length_hash, $tryfixrc);

around BUILDARGS => sub{
    my ($orig, $class, $line, $seqlengths, $tryfixrc) = @_;
    chomp($line);

    my ($readid, $flag, $seqid, $leftmost, $mapq, $cigar, $rnext, $pnext, $tlen, $readseq, $readqual, @optional) 
    = split /\t/, $line;

    my $mismatch_string;
    my $edit_distance = 0;
    my $mismatches_xm = 0;
    for my $o (@optional) {
        if ($o =~ /MD:Z:([A-Z\d\^]+)/){
            $mismatch_string = $1;
        }
        elsif ($o =~ /NM:i:(\d+)/){
            $edit_distance = $1;
        }
        elsif ($o =~ /XM:i:(\d+)/){
            $mismatches_xm = $1;
        }
    }

    croak "If \$tryfixrc is true, must also pass \$seqlengths" if $tryfixrc && ! $seqlengths;

    my $fixrc;
    if ($tryfixrc && $seqid =~ s/^RC_//){
        # if we get here, we know it's mapped since $seqid was not *
        $fixrc = 1;
        croak "\$tryfixrc is true, but \$seqid $seqid not found in \$seqlengths?"
        if (!exists $seqlengths->{uc $seqid});

        # $readid unchanged
        $flag = $flag & $flag_bits{reverse}  # is reverse?
              ? $flag & ~$flag_bits{reverse} # then unreverse it
              : $flag | $flag_bits{reverse}; # or reverse (there's probably a one-line way to do this...)

        # $seqid already stripped of RC_
        # $leftmost stored in original_leftmost and leftmost/rightmost/span are lazy.
        # $mapq unchanged
        # $cigar should be reversed, but should be done lazily since we probably won't need it
        # not sure about rnext 
        # not sure about pnext 
        # tlen probably unchanged 
        $readseq = reverse_complement($readseq);
        # $readqual = scalar reverse $readseq;
        $readqual = reverse_complement($readqual);
        
        # edit_distance can remain same
        # $mismatch_string should be reversed like cigar lazily.
    }

    return $class->$orig(
        readid   => $readid,
        flag     => $flag,
        seqid    => $seqid eq '*' ? undef : $seqid,
        original_leftmost => $leftmost,
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
        mismatches_xm            => $mismatches_xm,

        fixrc    => $fixrc,
        seqlen   => $tryfixrc && $seqlengths->{uc $seqid},
    );
};


#######################################################################
# accessors

# 11 mandatory fields are always there
has readid                => ( is => 'ro', required => 1 );
has flag                  => ( is => 'ro', required => 1 );
has seqid                 => ( is => 'ro', required => 1 ); #undef if unmapped
has original_leftmost     => ( is => 'ro', required => 1 );
has mapq                  => ( is => 'ro', required => 1 );
has original_cigar_string => ( is => 'ro', required => 1 );
has rnext                 => ( is => 'ro', required => 1 );
has pnext                 => ( is => 'ro', required => 1 );
has tlen                  => ( is => 'ro', required => 1 );
has readseq               => ( is => 'ro', required => 1 );
has readqual              => ( is => 'ro', required => 1 );

has original_mismatch_string => ( is => 'ro', required => 1,);
has edit_distance   => ( is => 'ro', required => 1,);
has mismatches_xm   => ( is => 'ro', required => 0,);

# true if tryfixrc and /^RC_/
has fixrc => ( is => 'ro', required => 1 );

# length of chromosome
has seqlen => ( is => 'ro', required => 1 );

sub readlength {    length($_[0]->readseq)                 }

# bit flag extraction
sub is_reverse {    $_[0]->flag & $flag_bits{reverse};     }
sub failed_qc  {    $_[0]->flag & $flag_bits{failed_qc};   }
sub mapped     { !( $_[0]->flag & $flag_bits{unmapped}   ) }
sub is_mapped  { !( $_[0]->flag & $flag_bits{unmapped}   ) }
sub is_first_segment       { $_[0]->flag & $flag_bits{first_segment} }
sub is_second_segment      { $_[0]->flag & $flag_bits{last_segment} }
sub is_primary_alignment   { !($_[0]->flag & $flag_bits{secondary_alignment}) }
sub is_secondary_alignment { $_[0]->flag & $flag_bits{secondary_alignment} }

# two alternatives for leftmost/rightmost when we know 
# the alignment has no gaps, like for bowtie1.
sub dumbstart { 
    my $self = shift;
    if ($self->fixrc){
        return $self->seqlen - ($self->original_leftmost + $self->readlength() - 1) + 1;
    }
    else{
        return $self->original_leftmost 
    }
} 
sub dumbend { 
    my $self = shift;
    if ($self->fixrc){
        return $self->seqlen - $self->original_leftmost + 1;
    }
    else{
        return $self->original_leftmost + $self->readlength() - 1
    }
} 

#######################################################################
# rightmost is the right most position in the genome that read aligns

# notes:
# rightmost = POS + #(M/D/=/X/N) - 1
# definitely not not s,h,p... 
# not I b/c that only advances read
# manual says N is intron (for mrna). 
# I think =/X are subtype of M?
#
# length of SEQ = M/I/S/=/X (from manual)

has rightmost => ( is => 'ro', lazy_build => 1 );
has leftmost  => ( is => 'ro', lazy_build => 1 );
has span => (is => 'ro', lazy_build => 1);

sub _build_span{
    my $self = shift;
    my $span;

    for my $c (@{$self->cigar}) {
        my ($type, $count) = @$c;
        if ($type ~~ [qw{M D = X N}]){
            $span += $count;
        }
    }
    return $span;
}

sub _build_leftmost{
    my $self = shift;
    return if (! $self->mapped);

    if ($self->fixrc){
        return $self->seqlen - ($self->original_leftmost + $self->span - 1) + 1;
    }
    else{
        return $self->original_leftmost;
    }
}

sub _build_rightmost{
    my $self = shift;
    return if (! $self->mapped);

    if ($self->fixrc){
        return $self->seqlen - ($self->original_leftmost) + 1;
    }
    else{
        return $self->original_leftmost + $self->span - 1;
    }
}

# for debugging.  supposed to be equal to read length
has cigarlength => ( is => 'ro', lazy_build => 1,);

sub _build_cigarlength{
    my $self = shift;
    my $cigar_length = 0;
    for my $c (@{$self->cigar}) {
        my ($type, $count) = @$c;

        if ($type ~~ [qw{M I = X S}]){
            $cigar_length += $count;
        }
    }
    return $cigar_length;

    #croak "BUG: length mismatch between what CIGAR tells us ($cigar_length) and length of SEQ (" . length($seq) . ")" 
    #if length($seq) != $cigar_length;
}

#######################################################################
# cigar string.  idiotic format.  returns:
# [ [ [MIDNSHP=X], COUNT ], ... ]

# manual = http://samtools.sourceforge.net/SAM1.pdf
# CIGAR: Description
# M alignment match (can be a sequence match or mismatch)
# I insertion to the reference (# bases were inserted with w.r.t. the genome)
# D deletion from the reference (# bases were deleted from the genome)
# N skipped region from the reference
# S soft clipping (clipped sequences present in SEQ)
# H hard clipping (clipped sequences NOT present in SEQ)
# P padding (silent deletion from padded reference)
# = sequence match
# X sequence mismatch

# 30M1D20M1I40M10N9M:           1234567890          1234567890*
# 1         11        21        31        41        51         61        71        81        91        101       111       121
# NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN*NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
#          GTTTTGTTGTTTATGAGGATGTTAAGATGGKTTTTGAGGTTAGCTTGTTGT AATATTAATAATTTTTTTTGGGTAACATTGAAAAAAATAGKKKKKKKKKKTTTTTTTTT
#          GTTTTGTTGTTTATGAGGATGTTAAGATGG*TTTTGAGGTTAGCTTGTTGTCAATATTAATAATTTTTTTTGGGTAACATTGAAAAAAATAG>>>>>>>>>>TTTTTTTTT
#                                        Deletion (1D)        Insertion (1I)                           Gap 9M

sub matched_chunks{
    my $self = shift;
    my $left = $self->leftmost;
    my $cigar = $self->cigar;
    my @accum;
    for my $c (@$cigar) {
        my ($type, $count) = @$c;
        if ($type eq 'M'){
            my $right = $left + $count - 1;
            push @accum, [$left, $right];
            $left = $right + 1;
        }
        elsif ($type eq 'I' || $type eq 'D' || $type eq 'N'){
            $left += $count;
        }
        else{
            croak "sorry, S/H/P/=/X unsupported currently";
        }
    }
    return @accum;
}

has cigar => ( is => 'ro', lazy_build => 1 );

sub _build_cigar{
    my $self = shift;
    return if ! $self->mapped;
    my $cigar = $self->original_cigar_string;

    my @accum;

    while ($cigar =~ /(\d+)([MIDNSHP=X])/g){
        my $count = $1;
        my $type = $2;
        push @accum, [$type, $count];
        # croak "I haven't added support for hard clipping yet, sorry"
        # if $type eq 'H';
    }
    if ($self->fixrc){
        # pretty sure this is wrong but g'nuff initially.
        @accum = reverse @accum;
    }
    return \@accum;
}

has cigar_string => (is => 'ro', lazy_build => 1);

sub _build_cigar_string{
    my $self = shift;
    # return '*' if ! $self->mapped;
    if (! defined $self->cigar){
        return '*';
    }
    if ($self->fixrc){
        # pretty sure this is wrong but g'nuff initially.
        return $self->original_cigar_string;
    }
    return join "", map { $_->[1] . $_->[0] } @{$self->cigar};
}

# hack to compensate for the improper handling to ref pos in MDZ string.
sub cigar_is_all_m{
    my $self = shift;
    return $self->original_cigar_string =~ /^\d+M$/;
}

#######################################################################
# mismatch string parseing (MD:Z: in opt fields) string.  
# (slightly less) idiotic format.  

# mismatch_tokens is the parse MD:Z: string, which no more information

# mismatches looks like:
# [ 
#      [ 'M', COUNT ] match for COUNT bases
#   or [ 'C', POSITION, BASE_IN_REF, BASE_IN_READ ] CHANGE. reference base is BASE, read base is something else
#   or [ 'D', POSITION, BASE_IN_REF ]               DELETION of BASE from reference.
# ]

# snps returns 
# [ 
#   [ POSITION, BASE_IN_REF, BASE_IN_READ ]  # just the 'C' from above.
# ]

# "The MD field aims to achieve SNP/indel calling without looking at the
# reference. For example, a string ‘10A5^AC6’ means from the leftmost reference
# base in the alignment, there are 10 matches followed by an A on the reference
# which is different from the aligned read base; the next 5 reference bases are
# matches followed by a 2bp deletion from the reference; the deleted sequence
# is AC; the last 6 bases are matches. The MD field ought to match the CIGAR
# string."

# 2013-09-02 - Wait, so when would cigar ever be useful? Answer: see bug note below

has mismatch_tokens => ( is => 'ro', lazy_build => 1 );

sub _build_mismatch_tokens { 
    my $self = shift; 
    my $fixrc = $self->fixrc; 

    my $mdstring = $self->original_mismatch_string;
    if (! defined $mdstring){
        return [];
    }

    my $in_deletion = 0;
    my $read = $self->readseq(); # already RC'd.

    my @accum;

    # M = match
    # C = change
    # D = deletion (from genome)

    while ($mdstring =~ m{
            (\d+) 
            | 
            (\^?)([A-Z])
        }xmg){
        # number means match of that many bases.
        if ($1){
            # $ref_position += $rcpos * $token;
            # $read_position += $rcpos * $token;

            push @accum, ['M', $1];
        }
        # carat means deletion
        elsif ($2){
            if ($fixrc){
                push @accum, ['D', reverse_complement $3];
            }
            else{
                push @accum, ['D', $3];
            }
            # $ref_position += $rcpos;
        }
        # no carat means change of base
        elsif ($3){
            if ($fixrc){
                push @accum, ['C', reverse_complement $3];
            }else {
                push @accum, ['C', $3];
            }
        }
    }

    if ($fixrc){
        # I think this is enough?
        @accum = reverse @accum;
    }
    return \@accum;
}

# 'C' only: [ [ POSITION, BASE_IN_REF, BASE_IN_READ ]  ... ]
sub snps{
    my $self = shift; 
    return [map { [@{$_}[1,2,3]] } grep { $_->[0] eq 'C' } @{$self->mismatches}]
}

# BUG:
# below does not work for gapped alignments-- only works when the cigar string is like
# 100M. the ref_position needs to be correlated with whether there have been gaps/clipping 
# in the CIGAR string.  Instead of trying to increment ref_position, need a subroutine
# read_pos_to_ref_pos which 'walks' the cigar string. i hate this format. 
has mismatches => ( is => 'ro', lazy_build => 1 );

sub _build_mismatches { 
    my $self = shift; 

    my $read = $self->readseq(); # already RC'd. 

    my $ref_position  = $self->leftmost; # see bugnote above
    my $read_position = 0;

    my @accum;

    # M = match
    # C = change
    # D = deletion (from genome)

    # mismatch_tokens also rc'd
    for my $token (@{$self->mismatch_tokens}) {
        # say STDERR "ref_position: $ref_position";
        # say STDERR "read_position: $read_position";
        my $type = $token->[0];
        if ($type eq 'M'){
            my $count = $token->[1];
            $ref_position += $count;
            $read_position += $count;
        }
        elsif ($type eq 'D'){
            my $deleted_bases = $token->[1];
            my $count = length($deleted_bases);
            $ref_position += $count;
        }
        elsif ($type eq 'C'){
            my $bases_in_reference = $token->[1];

            for my $base (split //, $bases_in_reference) {
                push @accum, ['C', $ref_position, $base, substr($read, $read_position, 1) ];
                $ref_position += 1;
                $read_position += 1;
            }
        }
        else{
            croak "impossible";
        }
    }
    # say STDERR "ref_position: $ref_position";
    # say STDERR "read_position: $read_position";

    return \@accum;
}

# rebuild MD:Z: field from tokens. if no fixrc, just the original.

has mismatch_string => ( is => 'ro', lazy_build => 1 );

sub _build_mismatch_string{
    my $self = shift;
    if ($self->fixrc){
        return join "", map 
        {
            my ($type, $token) = @$_;
            if ($type eq 'M'){
                $token;
            }
            elsif ($type eq 'D'){
                "^$token";

            }
            elsif ($type eq 'C'){
                $token;
            }
        } @{$self->mismatch_tokens};
    }
    else{
        return $self->original_mismatch_string;
    }
}

#######################################################################
# stringification

use overload '""' => \&stringify;

sub stringify{
    my $self = shift;

    return join("\t", 
        $self->readid,
        $self->flag,
        $self->seqid // '*', 
        $self->leftmost // '0',
        $self->mapq, 
        $self->cigar_string,
        $self->rnext, 
        $self->pnext, 
        $self->tlen, 
        $self->readseq, 
        $self->readqual, 
        ($self->mapped ?  ("MD:Z:" . $self->mismatch_string) : ()),
        ($self->mapped ?  ("NM:i:" . $self->edit_distance) : ()),
        # more options...
    );

# chr5:746182:746281:+:746192:746259	0	chr5	746182	255	100M	*	0	0	ATGAAAATAACATTTTTATATATATTATTTGTTGAAATAATTATAAATTTAGTATTTATATATGTATATATATTTTTCTTTATATTTTTTATATGTATAT	ATGAAAATAACATCTCCATATATATTATTTGCTGAAATAATCACAAACTCAGTATTCACATATGCATATACATCCCTCTCTATATTTCTTATATGTATAT	XA:i:2	MD:Z:10T66T22	NM:i:2
# chr4:178891:178990:-:178894:178948	0	RC_chr4	1400932	255	100M	*	0	0	ATATTGTTGATATTTTGTGAGATGTAATGATAAAAAATGAATCTTTATGTAATGATAAAAATAAAATATTATATTTTAAGTTTTGGGATTTTTATGCATT	ACATTGTCGATACCTTGTGAGATGCAATGATAAAAAACGAATCTTCATGTAATGATAAAAATAAAACACCATATTTCAAGTTTTGGGATCTTCATGCATT	XA:i:2	MD:Z:42T53T3	NM:i:2
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
