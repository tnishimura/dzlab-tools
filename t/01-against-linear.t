#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Test::More qw(no_plan);
use Tree::Range;
use POSIX;
use List::Util qw/min max/;

my $max = 10000;
my $N = 1000;

sub randcoord{ return int(rand($max))-$max/2; }

sub gen_rand_library{
    my @accum;
    for (1..$N){
        my $start = randcoord();
        my $end = randcoord();
        my $label = "$start -> $end";
        push @accum, [$start,$end, $label];
    }
    return \@accum;
}

sub search_linear{
    my ($library, $start, $end) = @_;
    my $query = [$start, $end];
    my @results;

    for my $range (@$library) {
        if (my $cover = _range_overlap($range, $query)){
            push @results, {
                overlap => $cover,
                item => $range->[2],
            };
        }
    }
    return @results;
}

sub tree_range_from_lib{
    my $lib = shift;

    my $sr = Tree::Range->new();

    for my $range (@$lib) {
        $sr->add(@$range);
    }
    $sr->finalize();
    return $sr;
}

sub sort_results{
    return [sort {$a->{item} cmp $b->{item}} @_];
}

# utility. return number of units overlapped but $x, $y. 
sub _range_overlap{
    my ($x,$y) = @_;
    my $start1 = min($x->[0], $x->[1]);
    my $end1   = max($x->[0], $x->[1]);

    my $start2 = min($y->[0], $y->[1]);
    my $end2   = max($y->[0], $y->[1]);

    if ($end1 >= $start2 && $end2 >= $start1){
        return min($end1, $end2) - max($start1, $start2)  + 1;
    }
    else {
        return 0;
    }
}

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

