package IntervalMerge;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use List::Util qw/max min/;
use List::MoreUtils qw/all uniq/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(interval_merge interval_deinterlace);

# input:
#   $list      - arrayref to list of elems to be merge
#   $get_start - sub to get start coord of $list elem
#   $get_end   - sub to get end coord of $list elem
# output:
#   [
#     [start1, end1, @list_of_elems_contributing_to_this_interval_1],
#     [start2, end2, @list_of_elems_contributing_to_this_interval_2],
#     ...
#   ]
sub interval_merge{
    my ($list, $get_start, $get_end) = @_;

    die "argument error" unless all { ref $_ eq 'CODE' } ($get_start, $get_end);
    die "argument error" unless ref $list eq 'ARRAY';

    # interval version of cmp/<=>
    #   obj1 is [start, end, @elems],
    #   obj2 is [start, end, $elem]
    # return:
    #   -1 if $obj1 is behind $obj2
    #   0  if $obj1 is overlapping $obj2
    #   1  if $obj1 is after $obj2
    my $cmp_by_interval = sub {
        my ($obj1, $obj2) = @_;
        my ($start1, $end1, $start2, $end2) = ($obj1->[0], $obj1->[1], $obj2->[0], $obj2->[1]);

        if ($end1 < $start2){
            return -1;
        }
        elsif ($end2 < $start1){
            return 1;
        }
        else{
            return 0;
        }
    };

    # merge interval, concat element list
    #   in_progess is [start, end, @elems]
    #   current is    [start, end, $elems]
    #   output is     [start, end, @elems]
    my $interval_merge = sub{
        my ($prev, $current) = @_;
        my ($start1, $end1, @accum) = @$prev;
        my ($start2, $end2, $newelem) = @$current;

        return [min($start1, $start2), max($end1, $end2), @accum, $newelem];
    };

    # convert lists to [start, end, element], sort by start
    my @sorted = sort { 
        $a->[0] <=> $b->[0]
    } 
    map { 
        [$get_start->($_), $get_end->($_), $_] 
    } @$list;

    my @merged;

    # in_progess is [start, end, @elems]
    my $in_progress = shift @sorted; # first one.

    for my $current (@sorted) {
        my $cmp = $cmp_by_interval->($in_progress, $current);
        if ($cmp == -1){ # if $current is complete past $in_progress
            push @merged, $in_progress;
            $in_progress = $current;
        }
        elsif ($cmp == 1){ # $current behind $in_progress shouldn't happen
            die "wtf? not sorted?";
        }
        else{ # overlapping-- merge to $in_progress
            $in_progress = $interval_merge->($in_progress, $current);
        }
    }
    push @merged, $in_progress; # last one
    return \@merged;
}

# turn:
# [              ]
#     [              ]
#         [   ]
# into:
# [  ][  ][   ][ ][  ]
#
# NOTE: $list cannot contain any gaps.  
#
# input:
#   $list      - arrayref to list of elems to be merge
#   $get_start - sub to get start coord of $list elem
#   $get_end   - sub to get end coord of $list elem
# output:
#   [
#     [start1, end1],
#     [start2, end2],
#     ...
#   ]
sub interval_deinterlace{
    my ($list, $get_start, $get_end) = @_;

    my %starts = map { $get_start->($_) => 1 } @$list;
    my %ends = map { $get_end->($_) => 1 } @$list;

    my @positions = uniq sort { $a <=> $b } keys(%starts), keys(%ends);
    my @deinterlaced;

    #say Dumper \@positions;

    for my $i (1 .. $#positions) {
        my $from = $positions[$i-1];
        my $to   = $positions[$i];

        my $from_is_start = exists $starts{$from};
        my $from_is_end   = exists $ends{$from};
        my $to_is_start   = exists $starts{$to};
        my $to_is_end     = exists $ends{$to};

        if (! $from_is_end && $from_is_start && $to_is_end){
            push @deinterlaced, [$from, $to]
        }
        elsif ($from_is_end && $from_is_start && $to_is_end){
            push @deinterlaced, [$from + 1, $to]
        }
        elsif ($from_is_start && $to_is_start){
            push @deinterlaced, [$from, $to - 1]
        }
        elsif ($from_is_end && $to_is_end){
            push @deinterlaced, [$from + 1, $to]
        }
        elsif ($from_is_end && $to_is_start && $to - $from > 1){
            push @deinterlaced, [$from + 1, $to - 1]
        }
        elsif ($from_is_end && $to_is_start && $from + 1== $to){
            # skip:
            # [   ][    ]
            #     ^^
        }
        else{
            croak "impossible case... bug?";
        }
    }

    return \@deinterlaced;
}

1;

