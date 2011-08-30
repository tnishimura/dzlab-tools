package GFF::Tree;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use Carp;
use autodie;    
use Tree::Range;
use GFF::Parser;

has tree => (
    traits    => ['Hash'],
    is        => 'ro',
    isa       => 'HashRef[Item]',
    default   => sub { {} },
    handles   => {
        set_tree     => 'set',
        get_tree     => 'get',
        has_no_trees => 'is_empty',
        has_tree     => 'exists',
        num_trees    => 'count',
        delete_tree  => 'delete',
        tree_keys    => 'keys',
        tree_values  => 'values',
    },
);

has file => (
    is => 'ro',
    required => 1,
);

has normalize => (
    is => 'ro',
    default => 1,
);

has lenient => (
    is => 'ro',
    default => 0,
);


sub BUILD{
    my ($self) = @_;
    my $p = GFF::Parser->new(file => $self->file, normalize => $self->normalize);
    while (defined(my $gff = $p->next())){
        my ($seq, $start, $end) = ($gff->sequence, $gff->start, $gff->end);

        if (! $self->has_tree($seq)){
            $self->set_tree($seq => Tree::Range->new());
        }
        $self->get_tree($seq)->add($start,$end,$gff);
    }
    for my $tree ($self->tree_values) {
        $tree->finalize();
    }
}

sub search_overlap{
    my ($self, $seq, $start, $end) = @_;
    if (! $self->has_tree($seq)){
        if (! $self->lenient){
            croak "no such sequence $seq in GFF::Tree";
        }
        else {
            return;
        }
    }

    return $self->get_tree($seq)->search_overlap($start,$end);
}

sub search{
    my ($self, $seq, $start, $end) = @_;
    if (! $self->has_tree($seq)){
        if (! $self->lenient){
            croak "no such sequence $seq in GFF::Tree";
        }
        else {
            return;
        }
    }

    return $self->get_tree($seq)->search($start,$end);
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;

