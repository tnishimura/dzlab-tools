package FastaReader;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use Carp;
use autodie;    

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
sub _get_length { $_[0]->length()->{uc $_[1]}; }
sub sequence_lengths { 
    my $self = shift;
    return map { $_ => $self->_get_length($_); } $self->sequence_list;
}

#######################################################################
# get or set entire sequence

has sequence => (
    traits    => ['Hash'],
    is        => 'ro',
    isa       => 'HashRef[Str]',
    default   => sub { {} },
);

sub _set_sequence { $_[0]->sequence()->{uc $_[1]} = $_[2]; }
sub _get_sequence { $_[0]->sequence()->{uc $_[1]}; }
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
sub _get_original_name { $_[0]->original_name()->{uc $_[1]}; }

sub sequence_list{
    my $self = shift;
    return sort values %{$self->original_name};
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
    default => sub { sub { s/^>// } },
    init_arg => 'ht',
    documentation => "sub which messes with header via \$_",
);

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

sub get_pretty{
    my ($self, $seqid, $start, $end, %opts) = @_;

    my $sequence = $self->get($seqid, $start, $end, %opts);

    my @accum = (">${seqid}_${start}_${end}");
    my $length = length $sequence;

    my $pos = 0;

    while ($pos < $length){
        push @accum, substr $sequence, $pos, 80; 
        $pos+=80;
    }
    return join "\n", @accum;
}

# coord = 'f' if coords rel to 5', 'r' if 3'
# base  = 1 or 0
# rc    = whether to rc chunk
sub get{
    my ($self, $seqid, $start, $end, %opt) = @_;
    my $coord     = defined $opt{coord} ? lc($opt{coord}) : 'f';
    my $rc        = $opt{rc} // ($coord eq 'r');
    my $base      = $opt{base} // 1;
    $seqid = uc $seqid;

    my $totlen = $self->_get_length($seqid);
    my $lastindex = $totlen - 1;

    if (! defined $start && ! defined $end ){
        if ($self->slurp()){
            my $whole = $self->_get_sequence($seqid);
            if ($rc){
                $whole =~ tr/acgtACGT/tgcaTGCA/;
                $whole = reverse $whole;
            }
            return $whole;
        }
        else{
            return $self->get($seqid, 1, $totlen);
            #croak "whole sequence get only supported with slurp()ing on";
        }
    }
    elsif (! defined $end){
        $end = $base + $lastindex;
    }
    elsif (! defined $start){
        $start = $base;
    }

    # everything in base 0 coord now.
    $start -= $base;
    $end   -= $base;

    if ($end < $start){
        croak "get: end ($end) < start ($start)?";
    }
    if ($end < 0 || $start < 0 || $end > $lastindex || $start > $lastindex ){
        croak "start/end = ($start/$end) out of bounds)";
    }
    if ($coord ne 'f' && $coord ne 'r'){
        croak "\$coord needs tobe 'f' or 'r', case insensitive ";
    }

    my $left;
    my $right;

    if ($coord eq 'r'){
        ($left,$right) = ($totlen - 1 - $end, $totlen - 1 - $start);
    }
    else {
        ($left,$right) = ($start,$end);
    }

    if ($self->slurp){
        my $full = $self->_get_sequence($seqid);
        my $retrieved = substr $full, $left, $right-$left +1;

        if ($rc){
            $retrieved =~ tr/acgtACGT/tgcaTGCA/;
            $retrieved = reverse $retrieved;
        }
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
        if ($rc){
            $sub =~ tr/acgtACGT/tgcaTGCA/;
            $sub = reverse $sub;
        }
        return $sub;
    }
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

=head2 sequence_lengths()

returns hash of sequence names to lengths.

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

=cut

