package Ends::NeighborMap;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Moose;
use Carp;
use autodie;    
use Tree::Range;
use GFF::Parser;
use BinarySearch qw/greatest_lower least_upper/;
use List::Util qw/min max/;
use overload '""' => \&stringify;
#######################################################################
# Upstream/Downstream finder

# parameters
has flag => (is => 'ro', default => 0,);
has distance => (is => 'ro', default => 5000,);
has prime => (is => 'ro', default => 5);
has flag_6_distance => (is => 'ro', default => 150);

has starts => ( is => 'rw',);
has ends   => ( is => 'rw',);

has interval => (
    traits    => ['Hash'],
    is        => 'ro',
    isa       => 'HashRef[Any]',
    default   => sub { {} },
    handles   => {
        set_interval  => 'set',
        get_interval  => 'get',
        num_intervals => 'count',
        all_id        => 'keys'
    },
);

has search_tree => (
    is => 'ro', 
    default => sub { Tree::Range->new() },
);

sub stringify{
    my $self = shift;
    my @accum;
    for my $id (sort $self->all_id) {
        my ($position, $strand, $overlapped, $flag_upstream, $flag_downstream)
        = $self->neighborhood($id);
        push @accum, join "\t", $id, $strand, $flag_upstream, $flag_downstream ;
    }
    return join "\n", @accum;
}

sub add{
    my ($self, $id, $strand, $start, $end) = @_;
    $self->set_interval($id, [$strand, $start, $end]);
}

sub finalize{
    my ($self) = @_;
    my @ints = values %{$self->interval};
    $self->starts([sort { $a <=> $b } map { $_->[1] } @ints]);
    $self->ends([sort { $a <=> $b } map { $_->[2] } @ints]);

    # using Tree::Range to look up if a particular end is overlapped by another
    # gene.  This is hackish b/c either this entire class should be part of
    # Tree::Range, use T::R as a primary data structure, or something similar,
    # instead of creating two data structures?  Would like to store all info in
    # T::R as the primary data structure, but T::R doesn't yet support unique
    # label for each item, so checking if position is start/end becomes hard... 
    # -TN, 2011-09-20
    my $tr = $self->search_tree();
    for my $int (@ints) {
        $tr->add($int->[1], $int->[2], 0);
    }
    $tr->finalize();
}

sub nearest_upstream_end{
    my ($self,$position) = @_;
    return greatest_lower($self->ends, $position);
}

sub nearest_downstream_start{
    my ($self, $position) = @_;
    return least_upper($self->starts, $position);
}

sub overlapped{
    my ($self, $position) = @_;
    # yes, this is on purpose, so I can compare with test deeply. 
    # don't know how to do that with general boolean values
    return scalar($self->search_tree()->search($position)) > 1 ? 1 : 0;
}

# Coordinates always w.r.t. 5' end of chromosome.
# my ($position, $strand, $overlapped, $flag_upstream, $flag_downstream)
# = $nm->neighborhood($id);
sub neighborhood{
    my ($self, $id) = @_;
    my ($strand, $start, $end) = @{$self->get_interval($id)};
    my $distance = $self->distance;
    my $prime = $self->prime;
    my $flag_distance = $self->flag_6_distance;
    my $length = $end - $start; # + 1 not necessary

    #say "$strand, $start ,$end, $distance, $prime, $flag_distance, $length";

    #######################################################################
    # Upstream end
    if ($strand eq '+' && $prime == 5 || $strand eq '-' && $prime == 3){
        my $overlapped = $self->overlapped($start);
        my $max_minus = $start - $distance;
        my $max_plus = $start + $distance;
        my $nearest_upstream_end = $self->nearest_upstream_end($start);
        my $nearest_downstream_start = $self->nearest_downstream_start($start);
        my $minus = $nearest_upstream_end // $max_minus;
        my $plus  = $nearest_downstream_start // $max_plus;

        #say "$max_minus $minus $plus $max_plus ($start, $end)";

        given ($self->flag){
            when (0){
                return ($start, $strand, $overlapped, $max_minus, $max_plus);
            }
            when (2){
                return ($start, $strand, $overlapped,
                    max($minus, $max_minus),
                    min($end, $plus, $max_plus),
                );
            }
            when (6){
                return ($start, $strand, $overlapped,
                    max($minus, $max_minus),
                    min($end - $flag_distance, $plus, $max_plus),
                );
            }
            when (7){
                return ($start, $strand, $overlapped,
                    defined $nearest_upstream_end ? 
                      max($nearest_upstream_end, $start - $length) : 
                      $start - $length,
                    defined $nearest_downstream_start ? 
                      min($end, $nearest_downstream_start, $start + $length) : 
                      min($end, $start + $length),
                    $length,
                );
            }
        }
    }
    #######################################################################
    # Downstream end
    else {
        my $overlapped = $self->overlapped($end);
        my $max_minus = $end - $distance;
        my $max_plus = $end + $distance;
        my $nearest_upstream_end = $self->nearest_upstream_end($end);
        my $nearest_downstream_start = $self->nearest_downstream_start($end);
        my $minus = $nearest_upstream_end // $max_minus;
        my $plus  = $nearest_downstream_start // $max_plus;

        given ($self->flag){
            when (0){
                return ($end, $strand, $overlapped, $max_minus, $max_plus);
            }
            when (2){
                return ($end, $strand, $overlapped,
                    max($start, $minus, $max_minus),
                    min($plus, $max_plus),
                );
            }
            when (6){
                return ($end, $strand, $overlapped,
                    max($start + $flag_distance, $minus, $max_minus),
                    min($plus, $max_plus),
                );
            }
            when (7){
                return ($end, $strand, $overlapped,
                    defined $nearest_upstream_end ? 
                        max($start, $nearest_upstream_end, $end-$length) : 
                        max($start, $end-$length), 
                    defined $nearest_downstream_start ? 
                        min($nearest_downstream_start, $end + $length) :
                        $end + $length,
                    $length,
                );
            }
        }
    }
    croak "shouldn't be here";
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;

=head1 SYNOPSIS

 Ends::NeighborMap->new(
     flag            => 0 | 2 | 6,
     distance        => 5000,
     prime           => 5 | 3,
     flag_6_distance => 150
 );

=head2 $map->add($id, $strand, $start, $end)

=head2 $map->finalize()

=head2 $map->neighborhood($id)

 my ($position, $strand, $overlapped, $flag_upstream, $flag_downstream)
 = $nm->neighborhood($id);

=cut
