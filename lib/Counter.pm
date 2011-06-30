package Counter;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use Carp;
use autodie;    

has '_counter' => (
    traits  => ['Counter'],
    is      => 'ro',
    isa     => 'Num',
    default => 0,
    handles => {
        _inc_counter   => 'inc',
        #dec_counter   => 'dec',
        #reset_counter => 'reset',
    },
);

has '_increment' => (
    is => 'ro',
    default => 10000,
    init_arg => 'increment',
);

has verbose => (
    is => 'ro',
    default => 0,
);

sub increment{
    my $self = shift;

    my $c = $self->_inc_counter;

    if ($c % $self->_increment == 0){
        say STDERR $c if $self->verbose;
        return $c;
    }
    return 0;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;

