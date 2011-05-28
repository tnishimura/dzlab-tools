package Tree::Range::Collection;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use Carp;
use autodie;    
use Tree::Range;

has trees => (
    is => 'ro',
    traits  => ['Hash'],
    isa     => 'HashRef[Tree::Range]',
    default => sub { { } },
    handles => {
        has_tree => 'exists',
        groups    => 'keys',
        set_tree => 'set',
        get_tree => 'get',
    },
);

sub add{
    my ($self,$tree_name,$start,$end,$object) = @_;
    if (! $self->has_tree($tree_name)){
        $self->set_tree($tree_name => Tree::Range->new());
    }
    $self->get_tree($tree_name)->add($start,$end,$object);
}

sub finalize{
    my ($self) = @_;
    for my $g ($self->groups) {
        $self->get_tree($g)->finalize();
    }
}
sub search_overlap{
    my ($self, $tree_name, $start, $end) = @_;

    if ($self->has_tree($tree_name)){
        return $self->get_tree($tree_name)->search_overlap($start,$end);
    }
    else{
        carp "search_overlap called on non-existant group $tree_name";
        return;
    }
}
sub search{
    my ($self, $tree_name, $start, $end) = @_;
    if ($self->has_tree($tree_name)){
        return $self->get_tree($tree_name)->search($start,$end);
    }
    else{
        carp "search called on non-existant group $tree_name";
        return;
    }
}
sub info{
    my $self = shift;
    my @accum;
    for my $g (sort $self->groups) {
        push @accum, sprintf("%s: %d leaves", $g, $self->get_tree($g)->num_leaves);
    }
    return join "\n",@accum;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;

