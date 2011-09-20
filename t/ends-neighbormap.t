#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More qw(no_plan);

BEGIN { use_ok( 'Ends::NeighborMap' ); }
require_ok( 'Ends::NeighborMap' );

my $id = "AAA";

{
    my $ap = Ends::NeighborMap->new();
    $ap->add(++$id,'+',10,20);
    $ap->add(++$id,'+',25,30);
    is($ap->num_intervals(), 2, "num intervals");

    $ap->add(++$id,'-',14,25);
    $ap->finalize();
    is($ap->nearest_upstream_end(30), 25, "upstream end");
    is($ap->nearest_downstream_start(14), 25, "downstream start");
}

SKIP: {
    #skip "not now";

    my $ap = Ends::NeighborMap->new();
    $ap->add(++$id,'+',6000,7000);
    $ap->add(++$id,'+',16000,17000);
    $ap->add(++$id,'+',1000,3000);
    $ap->finalize();
    is($ap->nearest_upstream_end(6000), 3000, "upstream end meow");
}

SKIP: {
    my $ap = Ends::NeighborMap->new(flag => 0);

    $ap->add('A','-',1000,3000);
    $ap->add('B','+',6000,7000);
    $ap->add('C','-',8000,10000);
    $ap->add('D','+',9000,11000);

    $ap->finalize();

    #is($ap->nearest_upstream_end(6000), 3000, 'upstream end');
    is_deeply([$ap->neighborhood('B')], [6000,'+', 0, 1000, 11000]);
    is_deeply([$ap->neighborhood('C')], [10000,'-',1, 5000, 15000]);
    #is_deeply([$ap->neighborhood(10000, 1)], [6000,0, 3000, 8000]);
}

