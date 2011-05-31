package FastaReader;
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

has filehandle => (
    is => 'rw',
);


has location => (
    traits    => ['Hash'],
    is        => 'ro',
    isa       => 'HashRef[Int]',
    default   => sub { {} },
    handles   => {
        set_location     => 'set',
        get_location     => 'get',
        has_location     => 'exists',
        list_locations   => 'keys'
    },
);

has length => (
    traits    => ['Hash'],
    is        => 'ro',
    isa       => 'HashRef[Int]',
    default   => sub { {} },
    handles   => {
        set_length   => 'set',
        get_length   => 'get',
        has_length   => 'exists',
        list_lengths => 'keys'
    },
);

has sequence => (
    traits    => ['Hash'],
    is        => 'ro',
    isa       => 'HashRef[Str]',
    default   => sub { {} },
    handles   => {
        set_sequence     => 'set',
        get_sequence     => 'get',
    },
);

has slurp => (
    is => 'ro',
    default => 0,
);

has normalize => (
    is => 'ro',
    default => 1,
);

has header_transform => (
    is => 'ro',
    default => sub { sub { s/^>// } },
    init_arg => 'ht',
    documentation => "sub which messes with header via \$_",
);

has start_position => (
    is => 'rw',
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
        croak "file argument to FastaReader needs to be file handle or file name" . Dumper $self;
    }

    $self->start_position(tell($fh));

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
            $self->set_location($current => tell $fh);
        }
        else{
            $lengths{$current}+=length $line;
            if ($self->slurp){
                push @{$sequences{$current}}, $line;
            }
        }
    }
    while (my ($seq,$len) = each %lengths) {
        $self->set_length($seq => $len);
        if ($self->slurp){
            $self->set_sequence($seq => join '', @{$sequences{$seq}});
        }
    }

    $self->filehandle($fh);
}

sub _get_iter{
    my ($self, $seq) = @_;
    if (! $self->has_location($seq)){
        croak "no such sequence $seq";
    }
    my $pos = $self->get_location($seq);
    my $fh = $self->filehandle();
    seek $fh, $pos, 0;

    my $done = 0;
    return sub{
        if (defined(my $line = <$fh>)){
            $line =~ tr/\r\n//d;
            if ($line =~ /^>/){
                $done = 1;
                return;
            }
            if ($self->normalize){
                return uc $line;
            }
            else{
                return $line;
            }
        }
        else{
            return;
        }
    };
}

# coord = 'f' if coords rel to 5', 'r' if 3'
# base  = 1 or 0
# rc    = whether to rc chunk
sub get{
    my ($self, $seqid, $start, $end, %opt) = @_;
    my $coord     = defined $opt{coord} ? lc($opt{coord}) : 'f';
    my $rc        = $opt{rc} // ($coord eq 'r');
    my $base      = $opt{base} // 1;
    $seqid = $self->normalize ? uc $seqid : $seqid;

    my $totlen = $self->get_length($seqid);
    my $lastindex = $totlen - 1;

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
        my $full = $self->get_sequence($seqid);
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

            else {
                croak "wha?";
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

