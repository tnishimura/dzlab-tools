package GFF::Tied;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Moose;
use Carp;
use autodie;    
use Tie::File;
use Fcntl 'O_RDONLY';
use GFF::Util;

has file => (
    is => 'ro',
    required => 1,
);

has tied_file => (
    is      => 'rw',
    init_arg => undef,
);

sub BUILD{
    my $self = shift;
    tie my @array, 'Tie::File', $self->file(), mode => O_RDONLY;
    $self->tied_file(\@array);
}

sub get{
    my ($self, $index) = @_;
    if ($index >= 0 && $index <= $#{$self->tied_file()}){
        parse_gff($self->tied_file()->[$index]);
    }
    else {
        return undef;
    }

}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
