#######################################################################
# Node Role

package Tree::Range::Node;
use strict;
use warnings;
use Moose::Role;
use List::Util qw/min max/;

# utility. return number of units overlapped but $x, $y. 
sub overlap{
    my ($self,$x,$y) = @_;
    my $start1 = $self->start;
    my $end1   = $self->end;

    my $start2 = min($x,$y);
    my $end2   = max($x,$y);

    if ($end1 >= $start2 && $end2 >= $start1){
        return min($end1, $end2) - max($start1, $start2)  + 1;
    }
    else {
        return 0;
    }
}

requires qw/start end midpoint is_leaf/;

#######################################################################
# Leaf Class

package Tree::Range::Leaf;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use Carp;
use autodie;    
use List::Util qw/min max/;

around BUILDARGS => sub{
    my ($orig, $class, $start, $end, $item) = @_;
    $start = min($start,$end);
    $end   = max($start,$end);
    return $class->$orig(start => $start, end => $end, item => $item, midpoint => ($start + $end)/2);
};

has start => ( is => 'ro', required => 1);
has end => ( is => 'ro', required => 1);
has item => ( is => 'ro', required => 1);
has midpoint => (is => 'ro', isa => 'Num', required => 1);

has is_leaf => (is => 'ro', default => 1);

sub to_string{
    my $self = shift;
    return sprintf("%s => %s (%s)", $self->start, $self->end, $self->item);
}

with 'Tree::Range::Node';

no Moose;
__PACKAGE__->meta->make_immutable;

1;

#######################################################################
# Internal Node Class

package Tree::Range::Internal;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use Carp;
use autodie;    
use List::Util qw/min max/;

has left => (
    is => 'ro',
    isa => 'Tree::Range::Node',
    required => 1,
);
has right => (
    is => 'ro',
    isa => 'Tree::Range::Node',
    required => 1,
);
has start => (
    is => 'ro',
    lazy_build => 1,
    isa => 'Num',
);
has end => (
    is => 'ro',
    lazy_build => 1,
    isa => 'Num',
);

has is_leaf => (is => 'ro', default => 0);

has midpoint => (is => 'ro', lazy_build => 1, isa => 'Num');

sub _build_start{
    my $self = shift;
    return min($self->left->start, $self->right->start);
}
sub _build_end{
    my $self = shift;
    return max($self->left->end, $self->right->end);
}

sub _build_midpoint{
    my $self = shift;
    return ($self->start + $self->end)/2;
}

with 'Tree::Range::Node';

no Moose;
__PACKAGE__->meta->make_immutable;

1;

#######################################################################
# Tree

package Tree::Range;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use Carp;
use autodie;    
use List::Util qw/min max/;

has leaves => (
    traits  => ['Array'],
    is      => 'ro',
    isa     => 'ArrayRef[Tree::Range::Node]',
    default => sub { [] },
    handles => {
        get_leaves => 'elements',
        add_leaf   => 'push',
        num_leaves => 'count',
    },
);

has root => (
    is => 'rw',
);

sub add{
    my ($self,$start,$end,$item) = @_;
    $self->add_leaf(Tree::Range::Leaf->new($start,$end, $item));
}

sub finalize{
    my $self = shift;

    # flip around the sorting every time so that when there is an odd number of items like 2**N + 1, 
    # that last item is included at a deeper level.
    my $odd = 1; 

    my @current_level = $self->get_leaves();
    
    while (@current_level > 1){
        my @next_level = ();
        @current_level = sort { $odd * ($a->midpoint <=> $b->midpoint) } @current_level;
        $odd *= -1;

        while (@current_level >= 2){
            push @next_level, Tree::Range::Internal->new(left => shift(@current_level), right => shift(@current_level));
        }
        if (@current_level){
            push @next_level, @current_level;
        }

        @current_level = @next_level;
    }

    $self->root($current_level[0]);
}

sub dump{
    my $self = shift;
    _dump($self->root);
}

sub _dump{
    my ($node,$level) = @_;
    $level //= 0;

    if (ref $node eq 'Tree::Range::Leaf'){
        printf("%s%d => %d [%s]\n", " " x $level, $node->start, $node->end, $node->item);
    }
    else{
        printf("%s%d => %d\n", " " x $level, $node->start, $node->end);
        _dump($node->left, $level+1);
        _dump($node->right, $level+1);
    }
}

sub search{
    my ($self, $start, $end) = @_;
    return map {$_->{item}} $self->search_overlap($start,$end);
}
sub search_overlap{
    my ($self, $start, $end) = @_;

    my @accum = ();
    _search_overlap($self->root, $start, $end, \@accum);

    return @accum;
}
sub _search_overlap{
    my ($node, $start, $end, $accum) = @_;

    if (my $o = $node->overlap($start,$end)){
        if ($node->is_leaf){
            push @$accum, {item => $node->to_string, overlap => $o};
        }
        else{
            _search_overlap($node->left, $start, $end, $accum);
            _search_overlap($node->right, $start, $end, $accum);
        }
    }
}

#sub _search_overlap{
#    my ($node, $start, $end) = @_;
#
#    if (my $o = $node->overlap($start,$end)){
#        if ($node->is_leaf){
#            say $node->to_string;
#        }
#        else{
#            _search_overlap($node->left, $start, $end);
#            _search_overlap($node->right, $start, $end);
#        }
#    }
#}

no Moose;
__PACKAGE__->meta->make_immutable;

1;





=head1 NAME

Tree::Range - The great new Tree::Range!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use Tree::Range;

    my $foo = Tree::Range->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 function1

=cut

sub function1 {
}

=head2 function2

=cut

sub function2 {
}

=head1 AUTHOR

<Author>, C<< <<email at localhost>> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-tree-range at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Tree-Range>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tree::Range

You can also look for information at:

=head1 LICENSE AND COPYRIGHT

Copyright 2011 <Author>.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.

=cut

1; # End of Tree::Range
