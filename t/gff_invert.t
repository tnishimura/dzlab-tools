#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw{
    any all none notall true false
    firstidx first_index lastidx last_index
    insert_after insert_after_string
    apply indexes
    after after_incl before before_incl
    firstval first_value lastval last_value
    each_array each_arrayref
    pairwise natatime
    mesh zip uniq distinct minmax part
};

use Test::More qw(no_plan);

require_ok("gff_invert.pl");

sub get_start{ return $_[0]->[0]; }
sub get_end{ return $_[0]->[1]; }

sub merge{
    return [min($_[0]->[0], $_[1]->[0]), max($_[0]->[1], $_[1]->[1])];
}

is_deeply(
    [ [0,15], [20,29], [30,40], ],
    bubble_merge(
        [ [5,15], [0,10], [20,29], [30,40], ],
        \&get_start, \&get_end, 
    ),
    "bubble_merge",
);

is_deeply(
    [ [0, 5], [15, 20], [29, 30], [40, 50]],
    invert(
        [ [5,15], [20,29], [30,40], ], 0, 50
    ),
    "invert",
);


