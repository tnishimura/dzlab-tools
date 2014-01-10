package CDZUtils::CountVec;
use v5.12.0;
use warnings FATAL => "all";
use Data::Dumper;
use Carp qw/carp croak confess cluck/;
use autodie;
use FFI::Raw;

# if there's a libcdzutils.so in the current dir, use it instead of having
# FFI::Raw search.  useful for testing.
my $libname = -f 'libcdzutils.so' ? "./libcdzutils.so" : 'libcdzutils.so';

# countvec_t * countvec_new(int size, int base, int initial);
my $countvec_new = FFI::Raw -> new(
  $libname, 'countvec_new',
  FFI::Raw::ptr, 
  FFI::Raw::int,
  FFI::Raw::int,
  FFI::Raw::double,
);

# void countvec_destroy(countvec_t ** v);
my $countvec_destroy = FFI::Raw -> new(
  $libname, 'countvec_destroy',
  FFI::Raw::void, 
  FFI::Raw::ptr, 
);

# double countvec_get(countvec_t * v, int pos);
my $countvec_get = FFI::Raw -> new(
  $libname, 'countvec_get',
  FFI::Raw::double, 
  FFI::Raw::ptr, 
  FFI::Raw::int, 
);

my $countvec_sizeof_double = FFI::Raw -> new( $libname, 'sizeof_double', FFI::Raw::int, );
my $countvec_sizeof_int    = FFI::Raw -> new( $libname, 'sizeof_int',    FFI::Raw::int, );

# int countvec_pack(countvec_t * v, void * vbuffer, int from, int to);
my $countvec_pack_double = FFI::Raw -> new(
  $libname, 'countvec_pack_double',
  FFI::Raw::int, 
  FFI::Raw::ptr, 
  FFI::Raw::ptr, 
  FFI::Raw::int, 
  FFI::Raw::int, 
);

# int countvec_pack(countvec_t * v, void * vbuffer, int from, int to);
my $countvec_pack_int = FFI::Raw -> new(
  $libname, 'countvec_pack_int',
  FFI::Raw::int, 
  FFI::Raw::ptr, 
  FFI::Raw::ptr, 
  FFI::Raw::int, 
  FFI::Raw::int, 
);

# int countvec_is_in_range(countvec_t * v, int from, int to);
my $countvec_is_in_range = FFI::Raw -> new(
  $libname, 'countvec_is_in_range',
  FFI::Raw::int, 
  FFI::Raw::ptr, 
  FFI::Raw::int, 
  FFI::Raw::int, 
);

# int countvec_increment_range(countvec_t * v, int from, int to, double value);
my $countvec_increment_range = FFI::Raw -> new(
  $libname, 'countvec_increment_range',
  FFI::Raw::int, 
  FFI::Raw::ptr, 
  FFI::Raw::int, 
  FFI::Raw::int, 
  FFI::Raw::double, 
);

# int countvec_multiply_range(countvec_t * v, int from, int to, double value);
my $countvec_multiply_range = FFI::Raw -> new(
  $libname, 'countvec_multiply_range',
  FFI::Raw::int, 
  FFI::Raw::ptr, 
  FFI::Raw::int, 
  FFI::Raw::int, 
  FFI::Raw::double, 
);

# int countvec_set_range(countvec_t * v, int from, int to, double value);
my $countvec_set_range = FFI::Raw -> new(
  $libname, 'countvec_set_range',
  FFI::Raw::int, 
  FFI::Raw::ptr, 
  FFI::Raw::int, 
  FFI::Raw::int, 
  FFI::Raw::double, 
);

my $countvec_first = FFI::Raw -> new(
  $libname, 'countvec_first',
  FFI::Raw::int, 
  FFI::Raw::ptr, 
);

my $countvec_last = FFI::Raw -> new(
  $libname, 'countvec_last',
  FFI::Raw::int, 
  FFI::Raw::ptr, 
);

#######################################################################

use Scalar::Util qw/looks_like_number/;

sub new {
    my $class =shift;
    my %opt = @_;
    my $size = delete $opt{size};
    my $base = delete $opt{base} // 0;
    my $initial = delete $opt{initial} // 0.0;
    my $self = bless {}, $class;

    if (! looks_like_number($size) || 
        ! looks_like_number($base) || 
        ! looks_like_number($initial) ||
        ! ($base == 0 || $base == 1)){
        croak __PACKAGE__ . " constructor bad args";
    }

    $self->{cv} = $countvec_new->call($size, $base, $initial);
    return $self;
}

sub in_range{
    my ($self, $from, $to) = @_;
    return $countvec_is_in_range->call($self->{cv}, $from, $to);
}

sub get_double {
    my ($self, $from, $to) = @_;
    $from //= $countvec_first->call($self->{cv});
    $to   //= $countvec_last->call($self->{cv});
    return if ! $self->in_range($from, $to);

    my $getsize = $to - $from + 1;
    my $num_bytes = $countvec_sizeof_double->call() * $getsize;

    my $mem_buf = FFI::Raw::MemPtr->new($num_bytes);
    $countvec_pack_double->call($self->{cv}, $mem_buf, $from, $to) or croak "asdf";
    my @vec = unpack "d" x $getsize, $mem_buf->tostr($num_bytes);
    return \@vec;
}

sub get_int {
    my ($self, $from, $to) = @_;
    $from //= $countvec_first->call($self->{cv});
    $to   //= $countvec_last->call($self->{cv});
    return if ! $self->in_range($from, $to);

    my $getsize = $to - $from + 1;
    my $num_bytes = $countvec_sizeof_int->call() * $getsize;

    my $mem_buf = FFI::Raw::MemPtr->new($num_bytes);
    $countvec_pack_int->call($self->{cv}, $mem_buf, 1, $getsize) or croak "asdf";
    my @accum = unpack "l" x $getsize, $mem_buf->tostr($num_bytes);
    return \@accum;
}

sub set_range {
    my ($self, $from, $to, $value) = @_;
    return if ! $self->in_range($from, $to);
    return $countvec_set_range->call($self->{cv}, $from, $to, $value);
}

sub multiply_range {
    my ($self, $from, $to, $value) = @_;
    return if ! $self->in_range($from, $to);
    return $countvec_multiply_range->call($self->{cv}, $from, $to, $value);
}

sub increment_range {
    my ($self, $from, $to, $value) = @_;
    return if ! $self->in_range($from, $to);
    return $countvec_increment_range->call($self->{cv}, $from, $to, $value);
}

1;

