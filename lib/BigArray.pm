package BigArray;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use PDL;
use Params::Validate qw/:all/;

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

    $self->{pdl}         = zeroes $type_dispatch{$opt{type}}, $opt{size};
    $self->{size}        = $opt{size};
    $self->{buffer_size} = $opt{buffer_size};
    $self->{base}        = $opt{base};
    $self->{buffer}      = [];

    return $self;
}

sub info{
    my $self = shift;
    return $self->{pdl}->info("sizeof(%T[%D]) =%M");
    #return $self->{pdl}->info("Type: %T Dim: %-15D State: %S, %M");
}

sub get{
    my ($self, $coord) = @_;
    $self->{pdl}->at($coord - $self->{base});
}

sub get_range{
    my ($self, $start, $end) = @_;
    $start -= $self->{base};
    $end -= $self->{base};
    $self->{pdl}->slice("$start:$end");
}

sub increment_range{
    my ($self, $start, $end, $val) = @_;
    $start -= $self->{base};
    $end -= $self->{base};
    $self->{pdl}->slice("$start:$end") += $val//1;
}

#######################################################################
# incrementing one base at a time

sub commit_increment{
    my $self = shift;
    my $indices = pdl $self->{buffer};
    indadd(1, $indices, $self->{pdl});
    undef(@{$self->{buffer}});
}

sub push_increment{
    my $self      = shift;
    my $buffer    = $self->{buffer};
    my $base      = $self->{base};

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
    if (@$buffer > $self->{buffer_size}){ 
        $self->commit_increment();
    }
}

1;

