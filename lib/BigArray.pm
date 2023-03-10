=head2 BigArray

 my $ba = BigArray(
                size        => NUM,
                type        => [ byte | short | ushort | long* | longlong | float | double ],
                base        => [ 0* | 1 ],
                buffer_size => 10000,
          );

increment given positions by 1 (buffered):

 $b->push_increment(1, 4, 3, 10); # these indices incremeted by one

increment by key/val pair (buffered):

 $b->push_pair(5, 123, 0, 454); # add 123 to position 5, 454 to pos 0.

increment a range by a value (not buffer):

 $b->increment_range($start, $end, $val);

get a position, a range, or entire thing (as a PDL).  All these commit buffers automatically:

 $b->get($position);
 $b->get_range($start, $end);
 $b->get_pdl();

commiting should be unnecessary.

=cut
package BigArray;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use PDL;
use Params::Validate qw/:all/;
use Hash::Util qw/lock_keys/;

$PDL::BIGPDL = 1;

our %type_dispatch = (
    byte     => byte(),
    short    => short(),
    ushort   => ushort(),
    long     => long(),
    longlong => longlong(),
    float    => float(),
    double   => double(),
);

sub new {
    my $class = shift;
    my $self = bless {}, $class;

    my %opt = validate(@_, {
            size => {
                type => SCALAR,
                regex => qr/^\d+$/,
                callbacks => {
                    'greater than zero' => sub { shift > 0 },
                },
            }, 
            buffer_size => {
                type => SCALAR,
                regex => qr/^\d+$/,
                default => 10000,
                callbacks => {
                    'greater than zero' => sub { shift > 0 },
                },
            }, 
            type => {
                type => SCALAR,
                regex => qr/^byte|short|ushort|long|longlong|float|double$/,
                default => 'long',
            },
            base => {
                type => SCALAR,
                regex => qr/^0|1$/,
                default => 0,
            },
        });

    $self->{pdl}         = zeroes $type_dispatch{$opt{type}}, $opt{size} + 1;
    $self->{size}        = $opt{size};
    $self->{buffer_size} = $opt{buffer_size};
    $self->{base}        = $opt{base};

    # buffer for push_increment/commit_increment
    $self->{buffer}      = [];

    # buffers for push_pair/commit_pair
    $self->{pbuffer_indices} = [];
    $self->{pbuffer_values}  = [];

    $self->{dirty_pair} = 0; 
    $self->{dirty_increment} = 0; 

    lock_keys %$self;

    return $self;
}

sub info{
    my $self = shift;
    return $self->{pdl}->info("sizeof(%T[%D]) =%M");
    #return $self->{pdl}->info("Type: %T Dim: %-15D State: %S, %M");
}

sub get{
    my ($self, $coord) = @_;
    $self->commit();
    $self->{pdl}->at($coord - $self->{base});
}

sub get_range{
    my ($self, $start, $end) = @_;
    $self->commit();

    $start -= $self->{base};
    $end -= $self->{base};
    $self->{pdl}->slice("$start:$end");
}

sub increment_range{
    my ($self, $start, $end, $val) = @_;
    $start -= $self->{base};
    $end -= $self->{base};
    # ??? - probably due to bug in sam parser, we get weird out of range errors
    if ($start >= $self->get_size || $end >= $self->get_size){
        warn "out of range: $start $end " . $self->get_size;
        return;
    }
    $self->{pdl}->slice("$start:$end") += $val//1;
}

sub get_pdl{ 
    my ($self) = @_;
    $self->commit;
    return $self->{pdl}; 
}

# 
sub get_size             { return $_[0]->{size}; }

# mostly for testing/debugging
sub get_buffer           { return $_[0]->{buffer}; }
sub get_buffer_residence { return scalar @{$_[0]->{ buffer}}; }
sub get_buffer_size      { return $_[0]->{buffer_size}; }

#######################################################################
# incrementing one base at a time
# 
# this is separate from *_pairs for two reasons:
# a) I didn't know indadd could handle pdls for values
# b) it's 30-40% faster, which is important for MethylCounter


sub commit_increment{
    my $self = shift;
    my $indices = pdl $self->{buffer};
    indadd(1, $indices, $self->{pdl});
    undef(@{$self->{buffer}});
    $self->{dirty_increment} = 0;
}

# return 1 if commited
sub push_increment{
    my $self   = shift;
    my $buffer = $self->{buffer};
    my $base   = $self->{base};

    for my $i (@_) {
        push @$buffer, $i - $base;

        # no array checking b/c pdl does it for us
        #if ($i >= $base && $i < $max_index){
        #    push @$buffer, $i + $base;
        #}
        #else{
        #    croak "increment out-of-bounds error!";
        #}
    }
    $self->{dirty_increment} = 1;

    if (@$buffer >= $self->{buffer_size}){ 
        $self->commit_increment();
        return 1; 
    }
    return 0;

}

#######################################################################
# commit_pair (the name "_pair" sucks... _pair maybe?) 
# 
# for adding arbitrary values at indexes pairs to pdl. 

sub commit_pair{
    my $self = shift;
    croak "commit_pair: vbufs not same size" . Dumper $self if @{$self->{pbuffer_values}} != @{$self->{pbuffer_indices}};
    my $values = pdl $self->{pbuffer_values};
    my $indices = pdl $self->{pbuffer_indices};
    #croak $values->dim(0) . " " . $indices->dim(0);
    indadd($values, $indices, $self->{pdl});
    undef(@{$self->{pbuffer_values}});
    undef(@{$self->{pbuffer_indices}});
    $self->{dirty_pair} = 0;
}

sub push_pair{
    my $self   = shift;
    my $base   = $self->{base};
    my $pbuffer_values = $self->{pbuffer_values};
    my $pbuffer_indices = $self->{pbuffer_indices};

    if (@_ % 2 != 0){
        croak "usage: \$ba->push_pairs ind, val, ind, val...";
    }

    while (@_){
        push @$pbuffer_indices, shift(@_) - $base;
        push @$pbuffer_values, shift(@_);
    }
    $self->{dirty_pair} = 1;

    if (@$pbuffer_indices >= $self->{buffer_size}){
        $self->commit_pair();
        return 1;
    }
    return 0;
}

#######################################################################

sub commit{
    my $self =shift;
    $self->commit_increment if ($self->{dirty_increment});
    $self->commit_pair if ($self->{dirty_pair});
}

#######################################################################
# create_collection

# create_collection->({
#   type => 'long',
#   base => 1,
#   buffer_size => 10000,
# },
# {
#   chr1 => 100,
#   chr2 => 300,
# }
sub create_collection{
    my ($constructor_options, $size_specs) = @_;

    return {
        map{
            my $size = $size_specs->{$_};
            $_ => BigArray->new(size => $size, %$constructor_options);
        } keys %$size_specs
    };
}


1;

