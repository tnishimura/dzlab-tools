#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Test::More tests => 10;
use Tree::Range;
use List::Util qw/min max/;
use FindBin;
use lib "$FindBin::Bin";
use testutils;



for (1..10){
    my @range = (randcoord,randcoord);

    my $lib = gen_rand_library();
    my $sr = tree_range_from_lib($lib);

    my $binary_results = sort_results($sr->search_overlap(@range));
    my $linear_results = sort_results(search_linear($lib, @range));

    #say Dumper $binary_results;
    #say Dumper $linear_results;

    if (!is_deeply($binary_results, $linear_results,"search_indices $_")){
        say "=========== Found a bug =============="; 


        say "===search range @range";
        say "===original library:";
        for my $r (@$lib) {
            my $start = min($r->[0], $r->[1]);
            my $end   = max($r->[0], $r->[1]);
            say "$start\t=>\t$end";
        }
        say "===\$sr->dump:";
        $sr->dump;
        say "===linear:";
        say Dumper $linear_results;
        say "===binary:";
        say Dumper $binary_results;
        die "ERROR";
    }
}

1;
