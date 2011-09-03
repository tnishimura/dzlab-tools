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

#######################################################################
# make_gff_tree

=head2 make_trees

  my ($numbins, $id_list_href, $trees_href) 
  = Ends::NeighborMap::make_trees(
        file => $file,
        prime => 3|5,
        flag => 0|2|6,
        distance => 5000,
        flag_6_distance => 5000,
        locus_tag => 'ID',
        binwidth => 100,
    );

=cut
sub make_trees{
    croak "uneven \%opt" if @_ % 2;
    my %opt = @_;
    my $file            = delete $opt{file} || croak "need file";
    my $prime           = delete $opt{prime} || croak "need prime";
    my $flag            = delete $opt{flag} // croak "need flag";
    my $distance        = delete $opt{distance} || croak "need distance";
    my $flag_6_distance = delete $opt{flag_6_distance} || croak "need flag_6_distance";
    my $locus_tag       = delete $opt{locus_tag} || croak "need locus_tag";
    my $binwidth        = delete $opt{binwidth} || croak "need binwidth";

    my $numbins = int(($distance * 2) / $binwidth);

    my %neighbors;
    my %id_list;
    my %trees;

    my $parser = GFF::Parser->new(file => $file);

    while (defined(my $gff = $parser->next())){
        my ($seq, $start, $end, $strand, $id) = ($gff->sequence, $gff->start, $gff->end, $gff->strand, $gff->get_column($locus_tag));
        #say $gff->to_string;
        next if ! defined $seq; 
        push @{$id_list{$seq}}, $id;
        if (! exists $neighbors{$seq}){
            $neighbors{$seq} =
                Ends::NeighborMap->new(
                    flag            => $flag,
                    distance        => $distance,
                    prime           => $prime,
                    flag_6_distance => $flag_6_distance,
                );
        }
        $neighbors{$seq}->add($id, $strand, $start, $end);
    }
    for my $seq (keys %id_list) {
        my $nm = $neighbors{$seq};
        my $tr = Tree::Range->new();
        $nm->finalize();

        for my $id (@{$id_list{$seq}}) {
            my ($position, $strand, $flag_upstream, $flag_downstream)
            = $nm->neighborhood($id);
            
            my $binnum = $strand eq '+' ? 0 : $numbins - 1; # 0->99 or 99->0
            my $start = $position - $distance + 1;
            my $end = $start + $binwidth - 1;
            # -299->-200, -199->-100, -99->0, 1->100, 101->200, 201->300
            
            while ($start <= $distance + $position){
                if (Tree::Range::_overlap($start, $end, $flag_upstream, $flag_downstream)){
                    $tr->add($start, $end, [$id, $binnum]); 
                }
                $binnum += $strand eq '+' ? 1 : -1;
                $start += $binwidth;
                $end += $binwidth;
            }
        }
        $tr->finalize();
        $trees{$seq} = $tr;
    }
    return ($numbins, \%id_list, \%trees);
}

1;
