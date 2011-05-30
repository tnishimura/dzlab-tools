#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use Tree::Range;
use Benchmark qw/timethese cmpthese/;

my $min = 1; 
my $max = 30429313;

sub rand_position{ return $min + int rand $max} 

my $tree = Tree::Range->new();

while (defined(my $line = <>)){
    chomp $line;
    my ($start, $end) = (split /\t/, $line)[3,4];
    $tree->add($start,$end,$line);
}

$tree->finalize();

cmpthese(-1,{ 
    linear => sub{
        $tree->_linear_search_overlap((rand_position) x 2);
    },
    binary => sub{
        $tree->search_overlap((rand_position) x 2);
    }
});

