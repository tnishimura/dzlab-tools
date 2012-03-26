#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use autodie;
use Test::More tests => 10;
use List::Util qw/min max/;
use Tree::Range;
use Tree::Range::TestUtils;

for (1..10){
    my @range = randrange();

    my $sr = gen_rand_tree();

    my $binary_results = sort_results($sr->search_overlap(@range));
    my $linear_results = sort_results($sr->_linear_search_overlap(@range));

    if (!is_deeply($binary_results, $linear_results,"search_indices $_")){
        print "=========== Found a bug ==============\n"; 

        print "===search range @range\n";
        print "===\$sr->dump:\n";
        $sr->dump;
        print "===linear:\n";
        print Dumper $linear_results;
        print "===binary:\n";
        print Dumper $binary_results;
        die "ERROR";
    }
}

1;
