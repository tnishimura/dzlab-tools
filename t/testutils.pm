package testutils;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;
use Tree::Range;
use List::Util qw/min max/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(_range_overlap sort_results tree_range_from_lib gen_rand_library search_linear randcoord);

my $max = 10000;
my $N = 5000;

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

1;
