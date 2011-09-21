#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use List::Util qw/sum/;

END {close STDOUT}
$| = 1;
use FindBin;
use lib "$FindBin::Bin/lib";

use Ends::NeighborMapCollection;

my $file = "t/test1.gff";
my $nmc = Ends::NeighborMapCollection->new(
    file            => $file,
    tag             => 'ID',
    flag            => 0,
    distance        => 5000,
    prime           => 5,
    flag_6_distance => 150,
    binwidth        => 100,
);
my $tab = build_table($nmc);
#add_to_table($tab, 'AT1G01910', 10, 100,100);
#add_to_table($tab, 'AT1G01910', 10, 110,200);
#say Dumper $tab;
#dump_table($tab);

# { id => [[0,0], [1,2], ... ] }
sub build_table{
    my ($nmc) = @_;
    return {
        map { 
            $_ => [map {[0,0]} (1 .. $nmc->numbins )] 
        } 
        $nmc->all_id
    };
}
sub add_to_table{
    my ($table, $id, $bin, $c, $t) = @_;
    $table->{$id}[$bin][0]+=$c;
    $table->{$id}[$bin][1]+=$t;
}
sub dump_table{
    my $table = shift;
    for my $id (sort keys %$table) {
        say join "\t", 
        $id, 
        map { 
            my ($c, $t) = @$_;
            $c + $t == 0 ? 'na' : $c/($c+$t);
        } @{$table->{$id}};
    }
}
