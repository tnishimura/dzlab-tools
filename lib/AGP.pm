package AGP;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use Carp;
use autodie;    
use List::Util qw/max/;


has file => (
    is => 'ro',
    required => 1,
);
has _groups => (
    is => 'rw',
    default => sub { { } },
);
has _lengths => (
    is => 'ro',
    default => sub { { } },
);
sub BUILD{
    my ($self) = @_;

    open my $fh, '<', $self->file;
    while (defined(my $line = <$fh>)){
        chomp $line;
        my ($group, $object_start, $object_end, $partnum, $type, $component_id, $component_start, $component_end, $orientation) =
        split /\t/, $line;
        croak unless defined $group && defined $object_start;

        push @{$self->_groups->{$group}}, {
            object_start     => $object_start,
            object_end       => $object_end,
            object_length    => $type ne 'N' ? $object_end - $object_start + 1 : undef,
            component_id     => uc $component_id,
            component_start  => $component_start,
            component_end    => $component_end,
            component_length => $type ne 'N' ? $component_end - $component_start + 1 : undef,
            orientation      => $orientation,
            type             => $type,
        };
    }
    close $fh;

    for my $group ($self->groups) {
        $self->_lengths->{$group} = max map { $_->{object_end} } @{$self->_groups->{$group}};
    }
}

sub groups{
    my $self = shift;
    return sort keys %{$self->_groups};
}

sub length{
    my $self = shift;
    my $group = shift;
    return $self->_lengths->{$group};
}

sub next{
    my $self = shift;
    my $group = shift;
    return pop @{$self->_groups->{$group}};
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;

