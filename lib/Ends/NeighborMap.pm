package Ends::NeighborMap;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Moose;
use Carp;
use autodie;    
use BinarySearch qw/greatest_lower least_upper/;
use List::Util qw/min max/;

#######################################################################
# Upstream/Downstream finder

has starts => ( is => 'rw',);
has ends   => ( is => 'rw',);
has flag => (is => 'ro', default => 0,);
has distance => (is => 'ro', default => 5000,);
has prime => (is => 'ro', default => 5);
has flag_6_distance => (is => 'ro', default => 150);

has interval => (
    traits    => ['Hash'],
    is        => 'ro',
    isa       => 'HashRef[Any]',
    default   => sub { {} },
    handles   => {
        set_interval     => 'set',
        get_interval     => 'get',
        #has_interval     => 'exists',
        num_intervals    => 'count',
        #delete_interval  => 'delete',
        #list_intervals   => 'keys'
    },
);

sub add{
    my ($self, $id, $strand, $start, $end) = @_;
    $self->set_interval($id, [$strand, $start, $end]);
}

sub finalize{
    my ($self) = @_;
    $self->starts([sort { $a <=> $b } map { $_->[1] } values %{$self->interval}]);
    $self->ends([sort { $a <=> $b } map { $_->[2] } values %{$self->interval}]);
}

sub nearest_upstream_end{
    my ($self,$position) = @_;
    return greatest_lower($self->ends, $position);
}

sub nearest_downstream_start{
    my ($self, $position) = @_;
    return least_upper($self->starts, $position);
}

=head2 position_info

Coordinates always w.r.t. 5' end of chromosome.

 my ($position, $strand, $flag_upstream, $flag_downstream)
 = $nm->neighborhood($id);

=cut
sub neighborhood{
    my ($self, $id) = @_;
    my ($strand, $start, $end) = @{$self->get_interval($id)};
    my $distance = $self->distance;
    my $prime = $self->prime;
    my $flag_distance = $self->flag_6_distance;

    #say "$strand, $start ,$end, $distance, $prime, $flag_distance";

    if ($strand eq '+' && $prime == 5 || $strand eq '-' && $prime == 3){
        my $max_minus = $start - $distance;
        my $max_plus = $start + $distance;
        my $minus = greatest_lower($self->ends, $start) // $max_minus;
        my $plus  = least_upper($self->starts, $start) // $max_plus;
        #say "$max_minus $minus $plus $max_plus ($start, $end)";


        given ($self->flag){
            when (0){
                return ($start, $strand, $max_minus, $max_plus);
            }
            when (2){
                return ($start, $strand, 
                    max($minus, $max_minus),
                    min($end, $plus, $max_plus),
                );
            }
            when (6){
                return ($start, $strand, 
                    max($minus, $max_minus),
                    min($end - $flag_distance, $plus, $max_plus),
                );
            }
        }
    }
    else {
        my $max_minus = $end - $distance;
        my $max_plus = $end + $distance;
        my $minus = greatest_lower($self->ends, $end) // $max_minus;
        my $plus  = least_upper($self->starts, $end) // $max_plus;

        given ($self->flag){
            when (0){
                return ($end, $strand, $max_minus, $max_plus);
            }
            when (2){
                return ($end, $strand, 
                    max($start, $minus, $max_minus),
                    min($plus, $max_plus),
                );
            }
            when (6){
                return ($end, $strand, 
                    max($start + $flag_distance, $minus, $max_minus),
                    min($plus, $max_plus),
                );
            }
        }
    }
    croak "shouldn't be here";
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
