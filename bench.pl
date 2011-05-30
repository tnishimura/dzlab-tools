#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use Tree::Range;

my $min = 1; 
my $max = 30429313;

sub rand_position{ return $min + int rand $max} 

my $tree = Tree::Range->new();
open my $fh, '<', 'TAIR8_genes.gff';

while (defined(my $line = <$fh>)){
    chomp $line;
    my ($start, $end) = (split /\t/, $line)[4,5];
    $tree->add($start,$end,$line);
}

close $fh;
$tree->finalize();

for (1 .. 1_000_000){
    $tree->search((rand_position) x 2);
}
