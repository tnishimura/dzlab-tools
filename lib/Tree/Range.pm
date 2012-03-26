#######################################################################
# Tree

package Tree::Range;
use strict;
use warnings;
use Data::Dumper;
use Moose;
use Carp;
use autodie;    
use Scalar::Util qw/looks_like_number/;
use 5.010;
use FindBin;
use lib "$ENV{HOME}/dzlab-tools/lib";


BEGIN {
    eval 'require Inline::C';

    my $pkg;
    if ($@) {
        require Tree::Range::Overlap_PP;
        $pkg = "Tree::Range::Overlap_PP";
    }
    else {
        require Tree::Range::Overlap_XS;
        $pkg = "Tree::Range::Overlap_XS";
    }
    no strict 'refs';
    my $class = __PACKAGE__;
    *{"$class\::_overlap"} = *{"$pkg\::_overlap"};
}

sub _min{ return ($_[0] < $_[1] ? $_[0] : $_[1]); }
sub _max{ return ($_[0] > $_[1] ? $_[0] : $_[1]); }

# Internals
# Leaf: [start, end, midpoint, 1, item]
# Internal Node: [start, end, midpoint, 0, left, right]
# The constructors should guarantee that start < end;

has leaves => (
    traits  => ['Array'],
    is      => 'ro',
    isa     => 'ArrayRef[Item]',
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

has finalized => (
    is => 'rw',
    default => 0,
);


sub _create_leaf{
    my ($start, $end, $item) = @_;

    if (!(defined $start && defined $end && defined $item && 
            looks_like_number($start) && looks_like_number($end))){
        croak "argument error in add()";
    }
    # round to integers
    $start = int($start+.5);
    $end = int($end+.5);

    ($start,$end) = (_min($start,$end), _max($start,$end));

    return [$start, $end, ($start + $end)/2, 1, $item];
}

sub _create_internal{
    my ($left, $right) = @_;

    my $start = _min($left->[0], $right->[0]);
    my $end   = _max($left->[1], $right->[1]);

    # no need to check since $left/$right were created with create_leaf and
    # were sanitized there (right?)

    return [$start, $end, ($start + $end)/2, 0, $left, $right];
}


sub add{
    my ($self,$start,$end,$item) = @_;
    if ($self->finalized()){
        croak "Cannot add after finalize()-ing";
    }
    $self->add_leaf(_create_leaf($start,$end, $item));
}

sub finalize{
    my $self = shift;
    if ($self->finalized()){
        carp "Finalized twice?";
    }

    # flip around the sorting every time so that when there is an odd number of items like 2**N + 1, 
    # that last item is included at a deeper level.
    my $odd = 1; 

    my @current_level = $self->get_leaves();
    
    while (@current_level > 1){
        my @next_level = ();
        @current_level = sort { $odd * ($a->[2] <=> $b->[2]) } @current_level;
        $odd *= -1;

        while (@current_level >= 2){
            push @next_level, _create_internal(shift(@current_level), shift(@current_level));
        }
        if (@current_level){
            push @next_level, @current_level;
        }

        @current_level = @next_level;
    }

    $self->finalized(1);
    $self->root($current_level[0]);
}

sub _dump{
    my $self = shift;
    _dump_helper($self->root);
}

sub _dump_helper{
    my ($node,$level) = @_;
    $level //= 0;

    if (ref $node eq 'ARRAY'){
        printf("%s%d => %d [%s]\n", " " x $level, $node->[0], $node->[1], $node->item);
    }
    else{
        printf("%s%d => %d\n", " " x $level, $node->[0], $node->[1]);
        _dump_helper($node->[4], $level+1);
        _dump_helper($node->[5], $level+1);
    }
}

sub search{
    my ($self, $start, $end) = @_;

    if (!$self->finalized()){
        croak "Cannot search until finalized";
    }

    if (!defined $start){
        croak "search needs at least a start";
    }
    $end //= $start; # allow querying for a single point

    return map {$_->{item}} $self->search_overlap($start,$end);
}

sub search_overlap{
    my ($self, $start, $end) = @_;

    if (!$self->finalized()){
        croak "Cannot search_overlap until finalized";
    }

    if (!defined $start && !defined $end){
        croak "search_overlap needs a range...";
    }
    # don't allow for single point here b/c there's no point... single point
    # will always be full overlap
    
    my @accum = ();
    _search_overlap_iter($self->root, $start, $end, \@accum);

    return @accum;
}

sub _search_overlap_iter{
    my ($node, $start, $end, $accum) = @_;

    my @search_queue = ($node);

    while (@search_queue){
        my $node = pop @search_queue;
        # Leaf: [start, end, midpoint, 1, item]
        # Internal Node: [start, end, midpoint, 0, left, right]
        my ($node_start, $node_end, $is_leaf, $left_or_item, $right) = @{$node}[0,1,3,4,5];

        #if (my $o = _overlap2($node->[0],$node->[1],$start,$end)){
        #if (my $o = _overlap($node->[0],$node->[1],$start,$end)){
        if (my $o = _overlap($node_start, $node_end,$start,$end)){
            if ($is_leaf){
                push @$accum, {item => $left_or_item, overlap => $o};
            }
            else{
                push @search_queue, $left_or_item, $right;
            }
        }
    }
}


# for debugging/benchmarking purposes
sub _linear_search_overlap{
    my ($self, $start, $end) = @_;

    my @accum;
    for my $l ($self->get_leaves()) {
        if (my $o = _overlap($l->[0],$l->[1],$start,$end)){
            push @accum, {item => $l->[4], overlap => $o};
        }
    }
    return @accum;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;

=head1 NAME

Tree::Range - A binary search tree for ranges.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

 use Tree::Range;

 my $tr = Tree::Range->new();
 
 $tr->add(10,20,"a");
 $tr->add(15,25,"b");
 $tr->add(50,67,"c");
 $tr->add(22,49,"d");

 $tr->finalize();

 $tr->search(19,23);

=head1 DESCRIPTION

Tree::Range implements a binary tree for ranges.  Given a list of ranges [x,y]
and associated items, it allows searching for all ranges which overlap with a
query range.  

The tree is constructed by sorting the leaf ranges by midpoints, joining
adjacent leaves into new nodes with a range just big enough to encompass both,
and repeating.  The idea is that each node fully covers any descendant nodes.

Search is done depth-first, only traversing nodes that overlap (at least
partially with the query range.

                                        |-----------------------|
 |-----------------|                   |----------|           |-------|
             |-----------------| |-------|                         

         |            |              |      |       |             |
         +------------+              +------+       +-------------+
               |                        |                  |
               |                        +------------------+
               |                                 |
               +---------------------------------+
                                |

Performance note: when there are no overlaps between the ranges, this data
structure provides O(log n) search time.  As the number of overlaps increase
and as the ranges clump closer together, performance 
degrades to O(N).  

Also note: this module currently assumes that the coordindates are integers,
and the width of a singleton interval [$n, $n] is 1, not 0.

=head1 SUBROUTINES/METHODS

=head2 add 

 $tr->add($start, $end, $item);

Add an $item on range [$start, $end]. $start and $end are rounded to the nearest integers.
Note that you must finalize() before searching.

=head2 finalize

 $tr->finalize()

Lock the tree from further add()-ing and build the tree internally. (TODO:
build tree incrementally instead).

=head2 search 

 $tr->search($start, $end)
 $tr->search($point)

Returns a list of $item's that you added with add(). Need to finalize() before calling.

=head2 search_overlap

 $tr->search_overlap($start, $end)

Returns a list of hashrefs, each with two keys: 'item', pointing to the
original objected added with add(), and 'overlap', which is the size of the
query overlapped by the result.  (For example, if [$start,$end] is [10,20] and
the tree has the range [19,21], overlap will be 2 (at 19 and 20)).  Need to
finalize() before calling.

 [{ item => $item, overlap => $size_of_query_overlapped }, ...]

=head1 DEPENDENCIES

=over 

=item *

Moose

=back

=head1 AUTHOR

Tom Burns, C<< <tom@burns.org> >>

=head1 BUGS

Probably...

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tree::Range

You can also look for information at:

=head1 LICENSE AND COPYRIGHT

Copyright 2011 <Author>.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

=cut

1; # End of Tree::Range


