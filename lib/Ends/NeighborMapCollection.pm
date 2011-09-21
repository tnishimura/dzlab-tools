package Ends::NeighborMapCollection;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Moose;
use Carp;
use autodie;    

has file            => ( is => 'ro', required => 1,);
has tag             => ( is => 'ro', required => 1, );
has flag            => (is  => 'ro', required => 1, );
has distance        => (is  => 'ro', required => 1, );
has prime           => (is  => 'ro', required => 1, );
has flag_6_distance => (is  => 'ro', required => 1, );
has binwidth        => (is  => 'ro', required => 1, );
has numbins         => (is  => 'rw', );

has id_list => (
    is        => 'rw',
    isa       => 'HashRef[Str]',
    default   => sub { {} },
);

has lookup_tree => (
    traits    => ['Hash'],
    is        => 'ro',
    default   => sub { {} },
    handles   => {
        set_lookup_tree     => 'set',
        get_lookup_tree     => 'get',
        has_lookup_tree     => 'exists',
    },
);

#sub add_score{
#    my ($hash, $id, $bin, $item) = @_;
#    push @{$hash{$id}{$bin}}, $item;
#}

sub BUILD{
    my ($self) = @_;

    my $binwidth = $self->binwidth;
    my $distance = $self->distance;
    my $flag = $self->flag;
    my $flag_6_distance = $self->flag_6_distance;
    my $prime = $self->prime;
    my $numbins = $self->numbins = int(($self->distance * 2) / $binwidth);

    my %neighbormaps;
    my %lookup_trees;
    my %id_list;


    my $parser = GFF::Parser->new(file => $self->file);

    while (defined(my $gff = $parser->next())){
        my ($seq, $start, $end, $strand, $id) = ($gff->sequence, $gff->start, $gff->end, $gff->strand, $gff->get_column($self->locus_tag));
        #say $gff->to_string;
        next if ! defined $seq; 
        push @{$id_list{$seq}}, $id;
        if (! exists $neighbormaps{$seq}){
            $neighbormaps{$seq} = 
            Ends::NeighborMap->new(
                flag            => $flag,
                distance        => $distance,
                prime           => $prime,
                flag_6_distance => $flag_6_distance,
            );
        }
        $neighbormaps{$seq}->add($id, $strand, $start, $end);
    }
    for my $seq (keys %id_list) {
        my $nm = $neighbormaps{$seq};
        my $tr = Tree::Range->new();
        $nm->finalize();

        IDLOOP:
        for my $id (@{$id_list{$seq}}) {
            my ($position, $strand, $overlapped, $flag_upstream, $flag_downstream)
            = $nm->neighborhood($id);

            next IDLOOP if $overlapped;
            
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
        $self->set_lookup_tree($seq,$tr);
    }
    $self->id_list(\%id_list);

}


no Moose;
__PACKAGE__->meta->make_immutable;

1;
