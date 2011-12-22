#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
#use List::Util qw(first max maxstr min minstr reduce shuffle sum);
#use List::MoreUtils qw{
#    any all none notall true false
#    firstidx first_index lastidx last_index
#    insert_after insert_after_string
#    apply indexes
#    after after_incl before before_incl
#    firstval first_value lastval last_value
#    each_array each_arrayref
#    pairwise natatime
#    mesh zip uniq distinct minmax part
#};

use Test::More qw(no_plan);

use IntervalMerge;

sub get_start{ return $_[0]->[0]; }
sub get_end{ return $_[0]->[1]; }

{
    my $elem1 = { start => 5, end  => 15, text => "elem1" };
    my $elem2 = { start => 0, end  => 10, text => "elem2" };
    my $elem3 = { start => 20, end => 29, text => "elem3" };
    my $elem4 = { start => 30, end => 40, text => "elem4" };
    is_deeply(
        [ 
        [0,15, $elem2, $elem1], 
        [20,29, $elem3], 
        [30,40, $elem4], 
        ],
        interval_merge(
            [ $elem1, $elem2, $elem3, $elem4 ],
            sub { $_[0]{start} },
            sub { $_[0]{end} },
        ),
        "interval_merge with elems",
    );
}

is_deeply(
    [ 
    [0,1], 
    [2,4], 
    [5,6], 
    ],
    interval_deinterlace(
        [ [0, 4], [2, 6], ],
        sub { $_[0][0] },
        sub { $_[0][1] },
    ),
    "interval_deinterlace 1",

);

is_deeply(
    [ 
    [0,4], 
    [5,10], 
    [11,15], 
    ],
    interval_deinterlace(
        [ 
        [0, 10],
        [5, 15],
        ],
        sub { $_[0][0] },
        sub { $_[0][1] },
    ),
    "interval_deinterlace 2",
);

is_deeply(
    [ 
    [0,4], 
    [5,10], 
    [11,15], 
    [16,20],
    ],
    interval_deinterlace(
        [ 
        [0, 10],
        [5, 15],
        [15, 20],
        ],
        sub { $_[0][0] },
        sub { $_[0][1] },
    ),
    "interval_deinterlace head-to-head",
);


