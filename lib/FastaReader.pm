package FastaReader;
use strict;
use warnings FATAL => "all";
use Data::Dumper;
use 5.010_000;
use Moose;
use Carp;
use autodie;    
use FindBin;
use lib "$FindBin::Bin/lib";
use DZUtil qw/reverse_complement c2t g2a/;

has filename => (
    is => 'ro',
    required => 1,
    init_arg => 'file',
);

has filehandle => (
    is => 'rw',
);

#######################################################################
# byte-position of start of sequence (not including header) of seq.

has location => (
    traits    => ['Hash'],
    is        => 'ro',
    isa       => 'HashRef[Int]',
    default   => sub { {} },
);

sub _set_location { $_[0]->location()->{uc $_[1]} = $_[2]; }
sub _get_location { $_[0]->location()->{uc $_[1]}; }

#######################################################################
# lengths of sequences

has length => (
    traits    => ['Hash'],
    is        => 'ro',
    isa       => 'HashRef[Int]',
    default   => sub { {} },
);

sub _set_length { $_[0]->length()->{uc $_[1]} = $_[2]; }
sub get_length { $_[0]->length()->{uc $_[1]}; }
sub sequence_lengths { 
    my $self = shift;
    my $mutator = shift;
    return map { 
        my $key = defined $mutator ? $mutator->($_) : $_;
        $key => $self->get_length($key); 
    } $self->sequence_list;
}

sub info{
    my ($self) = @_;
    my @ret;
    for my $seq (sort $self->sequence_list) {
        push @ret, "$seq\t" . $self->get_length($seq);
    }
    return join "\n", @ret;
}

#######################################################################
# get or set entire sequence

has sequence => (
    traits    => ['Hash'],
    is        => 'ro',
    isa       => 'HashRef[Str]',
    default   => sub { {} },
);

sub _set_sequence { $_[0]->sequence()->{uc $_[1]} = \$_[2]; }
sub _get_sequence { ${ $_[0]->sequence()->{uc $_[1]} } }
sub _get_sub_sequence {
    # always base 0, forward coord, rc = 0
    my ($self, $seqid, $left, $right) = @_;
    my $ref = $self->sequence()->{uc $seqid};
    return substr $$ref, $left, $right-$left +1;

    #my $retrieved = substr $self->_get_sequence($seqid), $left, $right-$left +1;
}
sub has_sequence { exists $_[0]->length()->{uc $_[1]}; }

#######################################################################
# Original names - all names are normalized to upper case, so need 
# a way to retrieve them if printing/etc

has original_name => (
    traits    => ['Hash'],
    is        => 'ro',
    isa       => 'HashRef[Str]',
    default   => sub { {} },
);

sub _set_original_name { $_[0]->original_name()->{uc $_[1]} = $_[2]; }
sub get_original_name { $_[0]->original_name()->{uc $_[1]}; }

sub sequence_list{
    my $self = shift;
    return sort values %{$self->original_name};
}

sub sequence_count{
    my $self = shift;
    my @seqs = $self->sequence_list;
    return scalar @seqs;
}

#######################################################################
# if true, entire fasta file read into memory and added to sequence() hash

has slurp => (
    is => 'ro',
    default => 0,
);

#######################################################################
# a regular expression which modifies $_ to get important part of header lines

has header_transform => (
    is => 'ro',
    # default => sub { sub { s/^>// } },
    # Get the first work after > only... hopefully this doesn't break anything
    default => sub { 
        sub { 
            my $orig = $_;
            s/^>//;
            my @parts = split ' ', $_;
            if (@parts == 0){
                croak "can't parse fasta header $orig";
            }
            else{
                $_ = $parts[0];
            }
        } 
    },
    init_arg => 'ht',
    documentation => "sub which messes with header via \$_",
);

#######################################################################
# Constructor

sub BUILD{
    my ($self) = @_;
    my $fh;
    if (defined $self->filename && -f $self->filename ){
        open $fh, '<', $self->filename
            or croak "cannot open $self->filename";
    }
    else {
        croak "file argument to FastaReader needs to be file handle or file name" . Dumper $self;
    }

    my $current;
    my %lengths; # use tmp hash b/c calling set_length every time is slow
    my %sequences;

    while (defined(my $line = <$fh>)){
        $line =~ tr/\r\n//d;
        if ($line =~ /^>(.*)$/){
            {
                local $_ = $line;
                $self->header_transform->();
                $current = $_;
            }
            my $original = $current;
            $current = uc $current;
            $self->_set_original_name($current,$original);

            $self->_set_location($current => tell $fh);
        }
        else{
            $lengths{$current}+=length $line;
            if ($self->slurp){
                push @{$sequences{$current}}, $line;
            }
        }
    }
    while (my ($seq,$len) = each %lengths) {
        $self->_set_length($seq => $len);
        if ($self->slurp){
            $self->_set_sequence($seq => join '', @{$sequences{$seq}});
        }
    }

    $self->filehandle($fh);
}

#######################################################################
# iterator

sub _get_iter{
    my ($self, $seq) = @_;
    $seq = uc $seq;
    if (! $self->has_sequence($seq)){
        croak "no such sequence $seq";
    }
    my $pos = $self->_get_location($seq);
    # may need to share fh in the future if there are thousands of seqs in
    # fasta and calling make_iterators()
    open my $fh, '<', $self->filename;
    seek $fh, $pos, 0; 

    #say "_get_iter on $seq from pos $pos $fh";

    my $done = 0;
    return sub{
        #seek $fh, $pos, 0; # seek in sub if sharing fh's.

        # this could be more sophisticated, like reading only X bytes at a time.
        if (defined(my $line = <$fh>)){ 
            $line =~ tr/\r\n//d;
            if ($line =~ /^>/){
                $done = 1;
                return;
            }
            # $pos = tell $fh; # necessary if seeking in sub{}
            return uc $line;
        }
        else{
            return;
        }
    };
}

sub make_iterators{
    my $self = shift;
    return map { $_ => $self->_get_iter($_) } $self->sequence_list();
}

#######################################################################
# get_pretty

sub get_pretty{
    my ($self, $name, $seqid, $start, $end, %opts) = @_;

    my $sequence = $self->get($seqid, $start, $end, %opts);

    my @accum = 
    defined $name ? ">$name" :
    defined $start && defined $end ? (">${seqid}_${start}_${end}") :
    defined $start ?  (">${seqid}_${start}_end") :
    defined $end ?  (">${seqid}_start_${end}") :
    (">${seqid}");

    my $length = length $sequence;

    my $pos = 0;

    while ($pos < $length){
        push @accum, substr $sequence, $pos, 80; 
        $pos+=80;
    }
    return join "\n", @accum, "\n";
}

#######################################################################
# get_context

=head2 get_context

allowable options:
  base          => 1
  undef_on_nonc => 0 # return undef on non c/g position (instead of croaking)
  dinucleotide  => 0 # return C{A,C,G,T,H}
  
for trinuc (default): 
  otherwise return context CG, CHG, CHH

for dinuc:
  return CH when triplet falls off edge
  otherwise, CG, CA, CT, CC, or CH (when not ACGT).

Automatically determine the strand, since this method only makes sense when 
first position is C.

=cut
sub get_context{
    my ($self, $seqid, $position, %opt) = @_;

    my $base          = delete $opt{base} // 1;
    my $undef_on_nonc = delete $opt{undef_on_nonc} // 0;
    my $dinucleotide  = delete $opt{dinuc} // 0;

    croak "get_context only takes options base, undef_on_nonc, dinucleotide" if keys %opt;

    my @split = split //, $self->get_context_raw( $seqid, $position, base => $base, rc => 0);

    my $reverse;
    if ($split[0] eq 'C'){
        $reverse = 0;
    }
    elsif ($split[0] eq 'G'){ 
        @split = split //, $self->get_context_raw( $seqid, $position, base => $base, rc => 1);
        $reverse = 1;
    }
    elsif ($undef_on_nonc){
        return;
    }
    else{
        croak "get_context called on non-C position: " 
        . join("", @split) . " (base = $base, pos = $position, seq = $seqid)";
    }

    if ($dinucleotide){
        if (@split < 2){ return 'CH'; }
        elsif ($split[1] eq 'A'){ return 'CA'; }
        elsif ($split[1] eq 'C'){ return 'CC'; }
        elsif ($split[1] eq 'G'){ return 'CG'; }
        elsif ($split[1] eq 'T'){ return 'CT'; }
        else{ return 'CH'; }
    }
    else{
        if (@split >= 2 && $split[1] eq 'G'){ return 'CG'; }
        elsif (@split == 3 && $split[2] eq 'G'){ return 'CHG'; }
        elsif (@split == 3){ return 'CHH'; }
        else { return 'CHH'; }
    }
}

# return (up to) 3 bases from position. Not sure if this has utility
# except for get_context.
sub get_context_raw{
    my ($self, $seqid, $position, %opt) = @_;

    # extract options
    $seqid = uc $seqid;
    my $rc        = $opt{rc} // 0;
    my $base      = $opt{base} // 1;
    $opt{lenient} = 1;

    my ($start, $end) = $rc ? ($position - 2, $position) : ($position, $position + 2);
    return uc $self->get($seqid, $start, $end, %opt);
}
#######################################################################
# get - main sequence retrieval function

# coord = 'f' if coords rel to 5', 'r' if 3'
# base  = 1 or 0
# rc    = whether to rc chunk
# bs    = c2t, g2a
sub get{
    my ($self, $seqid, $start, $end, %opt) = @_;
    my $coord     = defined $opt{coord} ? lc($opt{coord}) : 'f';
    my $rc        = $opt{rc} // ($coord eq 'r');
    my $bs        = $opt{bs};
    my $base      = $opt{base} // 1;
    $seqid = uc $seqid;

    if (! $self->has_sequence($seqid)){
        croak "no such sequence $seqid";
    }

    if (defined $bs && $bs ne 'c2t' && $bs ne 'g2a'){
        croak "bs, if given, should be c2t or g2a";
    }

    my $totlen = $self->get_length($seqid);
    my $lastindex = $totlen - 1;

    if (! defined $start && ! defined $end ){
        my $whole = $self->slurp() ? $self->_get_sequence($seqid) : $self->get($seqid, 1, $totlen);

        $whole = reverse_complement($whole) if $rc;
        $whole = c2t $whole if defined $bs && $bs eq 'c2t';
        $whole = g2a $whole if defined $bs && $bs eq 'g2a';
        return $whole;
    }
    elsif (! defined $end){
        $end = $base + $lastindex;
    }
    elsif (! defined $start){
        $start = $base;
    }

    ### 
    # everything in base 0 coord now.
    $start -= $base;
    $end   -= $base;

    if ($opt{lenient}){
        if ($start < 0){ $start = 0; }
        if ($lastindex < $end){ $end = $lastindex; }
    }

    ###
    # error check
    if ($end < $start){
        croak "get: end ($end) < start ($start)?";
    }
    if ($end < 0 || $start < 0 || $end > $lastindex || $start > $lastindex ){
        croak "start/end = ($start/$end) out of bounds)";
    }
    if ($coord ne 'f' && $coord ne 'r'){
        croak "\$coord needs to be 'f' or 'r', case insensitive ";
    }

    ###
    # left/right are the 0-base, forward strand coordinates.
    my $left;
    my $right;

    if ($coord eq 'r'){
        ($left,$right) = ($totlen - 1 - $end, $totlen - 1 - $start);
    }
    else {
        ($left,$right) = ($start,$end);
    }

    if ($self->slurp){
        my $retrieved = $self->_get_sub_sequence($seqid, $left, $right);
        #my $retrieved = substr $self->_get_sequence($seqid), $left, $right-$left +1;
        #warn $retrieved;

        $retrieved = reverse_complement($retrieved) if $rc;
        $retrieved = c2t $retrieved if defined $bs && $bs eq 'c2t';
        $retrieved = g2a $retrieved if defined $bs && $bs eq 'g2a';
        return $retrieved;
    }
    else{

        #say STDERR "=== ($left, $right)";

        my $iter = $self->_get_iter($seqid);

        my $chunk_first_pos = 0;
        my @accum;

        while (defined(my $chunk = $iter->())){
            my $chunk_length = length $chunk;
            my $chunk_last_pos =  $chunk_first_pos + $chunk_length - 1;

            #                          left            right
            #                          |---------------|
            # |---------------| chunk
            if ($chunk_last_pos < $left){
                #warn "after $chunk_first_pos, $chunk_last_pos ($left, $right)";
            }

            # left            right
            # |---------------|
            #                     |---------------| chunk
            elsif ($right < $chunk_first_pos){
                #warn "before $chunk_first_pos, $chunk_last_pos ($left, $right)";
                last;
            }

            #           left            right
            #           |---------------|
            # |---------------| chunk
            elsif ($chunk_first_pos <= $left && $chunk_last_pos < $right ){
                #warn "overlap end of chunk $chunk_first_pos, $chunk_last_pos ($left, $right)";
                push @accum, substr $chunk, $left - $chunk_first_pos;
            }

            # left            right
            # |---------------|
            #          |---------------| chunk
            elsif ($left < $chunk_first_pos && $right <= $chunk_last_pos){
                #warn "overlap start of chunk $chunk_first_pos, $chunk_last_pos ($left, $right)";
                push @accum, substr $chunk, 0, $right - $chunk_first_pos + 1;
            }

            #     left  right
            #     |-----|
            # |---------------| chunk
            elsif ($chunk_first_pos <= $left && $right <= $chunk_last_pos){
                #warn "full overlap $chunk_first_pos, $chunk_last_pos ($left, $right)";
                push @accum, substr $chunk, $left - $chunk_first_pos, $right-$left +1;
            }

            # left              right
            # |-----------------|
            #     |----------|    chunk
            elsif ($left <= $chunk_first_pos && $chunk_last_pos <= $right){
                push @accum, $chunk;
            }

            else {
                croak "wha? chunk: ($chunk_first_pos, $chunk_last_pos), left/right: ($left, $right)";
            }

            $chunk_first_pos = $chunk_first_pos + $chunk_length;
        }

        my $sub = join '', @accum;

        $sub = reverse_complement($sub) if $rc;
        $sub = c2t $sub if defined $bs && $bs eq 'c2t';
        $sub = g2a $sub if defined $bs && $bs eq 'g2a';

        return $sub;
    }
}

#######################################################################
# Utilities

# given a sequence and 5' coordinate, return the 3' coordinate 
sub reverse2forward{
    my ($self, $seqid, $coord, $base) = @_;
    $base //= 1;

    $coord -=$base; # zero base
    my $totlen = $self->get_length($seqid);

    croak "$seqid does not exist" if ! $self->has_sequence($seqid);
    croak "$coord out of range" if ($coord < 0 && $totlen <= $coord);

    return $totlen - 1 - $coord + $base;
}

# they are actually the same... aren't they?
sub forward2reverse{ 
    return reverse2forward(@_);
}

=head2 format_fasta('header', $seq)

format a header and a single seq for printing as a fasta.
(Doesn't belong in FastaREADER module... but oh well)

=cut 
sub format_fasta{
    my ($header, $seq, $width) = @_;

    confess "need a sequence and header" unless ($seq and $header);

    my @buffer;
    $buffer[0] = ">$header" if defined $header;
    $width //= 80;
    my $length = length $seq;

    for (my $position = 0; $position < $length; $position += $width){
        push @buffer, substr $seq, $position, $width;
    }
    return (join "\n", @buffer) . "\n";
}

our @iupac_bases = qw/A C G T R Y S W K M B D H V N/;

our $iupac_bases_regex = qr/^[ACGTRYSWKMBDHVN]+$/;

# deprecated by Fasta::BaseComposition, but keep for testing purposes
sub base_composition{
    my $self = shift;

    # this is a hash so we can check existance fast:
    my %alpha = map { $_ => 0 } @iupac_bases, 'other';

    my $bufsize = 1_000_000;
    my %seqlens = $self->sequence_lengths;

    my %bases = map { $_ => { map { $_ => 0 } keys %alpha } }  keys %seqlens;

    while (my ($seq,$len) = each %seqlens) {
        my $position = 0;
        while ($position < $len) {
            #say "$seq, $len,  $position";
            my $subsec = $self->get(
                $seq, 
                $position, 
                $position + $bufsize -1, 
                lenient => 1, 
                base => 0,
            ); 
            for my $base (split //, $subsec){
                $base = uc $base;
                if (exists $alpha{$base}){
                    $bases{$seq}{$base}++;
                }
                else{
                    $bases{$seq}{other}++;
                }
            }
            $position+= $bufsize;
        }
    }
    return \%bases;
}

sub base_composition_table{
    my $self = shift;

    my $base_compo = $self->base_composition();

    my %seqlens = $self->sequence_lengths;

    my @accum;
    push @accum, join "\t", qw/Seq Length/, @iupac_bases, 'other';

    for my $seq (sort keys %seqlens) {
        push @accum, join "\t", $seq, $seqlens{$seq}, @{$base_compo->{$seq}}{@iupac_bases, 'other'};
    }
    return join "\n", @accum;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;



=head1 NAME
 
FastaReader - Fasta reader
 
=head1 VERSION
 
This documentation refers to FastaReader version 0.0.1
 
=head1 SYNOPSIS
 
    use FastaReader;

    my $f = FastaReader->new(file => 'file.fasta', ht => sub { s/>(\w+)/$1/; return $_ }, slurp => 0 );
    say $f->get('chr1', 1, 100);
  
=head1 DESCRIPTION
 
Lazy Fasta reader.

=head1 SUBROUTINES/METHODS 

=head2 FastaReader->new(file => 'filename or handle', slurp => 0 | 1, ht => sub {s/>(\w+)/$1/; return $1; })

=head2 get_length(SEQNAME)

returns length of specific sequence. Use this over sequence_lengths() hash, since this will
take care of name normalization.

=head2 sequence_list()

returns list of sequence names, sorted.

=head2 has_sequence()

check if particular sequence was in fasta file.

=head2 get(SEQNAME, START, END, coord => 'f' | 'r', base => 0 | 1, rc => 0 | 1)

Get subsequence (START, END) from SEQNAME.  coord is from forward strand by
default, but also be in reverse strand coordinates.  'rc => 1' reverse
compliments the retrieved subsequence. base is the coordinate of the first
base.

Default is coord => 'f', base => 1, rc => 0 (if coord => 'f') or 1 (if coord =>
'r').

=head2 get_pretty(...)

Exactly the same as get() except it returns a FASTA-formated string (with
header line and line breaks), suitable for printing.

=head2 get_context(SEQNAME, POS, base => 0 | 1, rc => 0 | 1)

Get the context at position POS.  Return 'CG', 'CHH', or 'CHG', or undef when
POS is an edge case and falls off edge without proper context. Croak when POS
out of bounds.

=cut

