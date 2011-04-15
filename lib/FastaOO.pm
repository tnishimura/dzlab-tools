package FastaOO;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use Carp;
use autodie;    

has filename_or_handle => (
    is => 'ro',
    required => 1,
    init_arg => 'file',
);

has normalize => (
    is => 'ro',
    default => 1,
);

has fasta => (
    traits    => ['Hash'],
    is => 'rw',
    isa       => 'HashRef[Str]',
    handles   => {
        get_full  => 'get',
        set_full  => 'set',
        seqs => 'keys',
    },
    default => sub { {} },
);

sub BUILD{
    my ($self) = @_;
    my $fh;
    if (ref $self->filename_or_handle eq 'GLOB'){
        $fh = $self->filename_or_handle;
    }
    elsif (!ref $self->filename_or_handle && -f $self->filename_or_handle ){
        open $fh, '<', $self->filename_or_handle
            or croak "cannot open $self->filename_or_handle";
    } elsif (! -f $self->filename_or_handle){
        croak $self->filename_or_handle . " doesn't exist?";
    }
    else {
        croak "file argument to GFF::Parser needs to be file handle or file name" . Dumper $self;
    }
    $self->slurp($fh);
    if (!ref $self->filename_or_handle && -f $self->filename_or_handle ){
        close $fh;
    }
    
}

=head1 EXPORTED FUNCTIONS

=head2 slurp_fasta('text.fasta', normalize => 1);

return hash of seqid (the first word after the '>') to the sequence. 
All sequence names upper cased if normalize => 1 (default)

=cut

sub slurp {
    my ($self,$fh) = @_;

    my %accum = ();

    my $current;
    my @buffer;
    while (defined(my $line = scalar readline $fh)) {
        $line =~ tr/\r\n//d;
        if ($line =~ /^>(\w+)/){
            if ($current){
                $self->set_full($current => join(q{}, @buffer));
                @buffer = ();
            }
            $current = $self->normalize() ? uc $1 : $1;
        } 
        else{
            push @buffer, $self->normalize() ? uc $line : $line;
        }
    }

    # last one
    if ($current){
        my $last = join(q{}, @buffer);
        $last =~ s/\s*//g;
        $self->set_full($current => $last);
    }
}

=head2 $self->get($seqid, $start, $end, coord => 'f', rc => 1, base => 0, normalize => 1)

If coord => 'r', use coordinates with respect to the 3' end. 
If rc => 1 (default for reverse coord), reverse complement. 

=cut

sub get{
    my ($self, $seqid, $start, $end, %opt) = @_;
    my $coord     = defined $opt{coord} ? lc($opt{coord}) : 'f';
    my $rc        = $opt{rc} // ($coord eq 'r');
    my $base      = $opt{base} // 1;
    my $normalize = $self->normalize;
    $seqid = $normalize ? uc $seqid : $seqid;

    if ($end < $start){
        croak "fasta_get: end ($end) < start ($start)?";
    }
    if ($coord ne 'f' && $coord ne 'r'){
        croak "FastaOO::Get: coord needs tobe 'f' or 'r', case insensitive ";
    }

    # everything in base 0 coord now.
    $start -= $base;
    $end   -= $base;
    my $sublen = ($end-$start)+1;

    if(my $seq = $self->get_full($seqid)){
        my $totlen = length $seq;
        my $lastindex = $totlen - 1;
        my $sub;
        if ($coord eq 'f'){
            my $left = $start;
            my $right = $left + $sublen - 1; 

            croak "($start,$end) (left = $left, right = $right) out of bounds" if ($left < 0 || $right > $lastindex);
            $sub = substr($seq,$left, $sublen);
        } 
        elsif ($coord eq 'r') {
            my $left = $totlen - $end - 1; 
            my $right = $left + $sublen - 1;
            croak "($start,$end) (left = $left, right = $right) on reverse out of bounds" if ($left < 0 || $right > $lastindex);

            $sub = substr($seq,$left, $sublen);
        }
        else {
            croak "coord needs to be 'f' or 'r'";
        }
        if ($rc){
            $sub =~ tr/acgtACGT/tgcaTGCA/;
            $sub = reverse $sub;
        }
        return $sub;
    } else {
        croak "no such seqid $seqid";
    }
}

=head2 $self->count()

Return hash of {seqid => count}.

=cut
sub counts {
    my ($self) = @_;
    return {map { $_ => length $self->get_full($_) } $self->seqs};
}

sub count{
    my ($self,$seqid) = @_;
    $seqid = $self->normalize ? uc $seqid : $seqid;
    return length $self->get_full($seqid);
}


=head2 format_fasta('header', $seq)

format a header and a single seq for printing as a fasta.

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

=head2 write($seqs)

Print fasta to file or file handle.

=cut

sub write{
    my ($self, $filename_or_handle) = @_;
    my $fh;
    if (ref $filename_or_handle eq 'GLOB'){
        $fh = $filename_or_handle;
    } else{
        open $fh, '>', $filename_or_handle;
    }

    for my $seq (sort $self->seqs) {
        say $fh format_fasta(uc $seq, $self->get_full($seq));
    }
    if (ref $filename_or_handle ne 'GLOB'){
        close $fh;
    }
}


no Moose;
__PACKAGE__->meta->make_immutable;

1;
