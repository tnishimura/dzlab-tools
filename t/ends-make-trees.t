package ends::make::trees;
use base qw(Test::Class);
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More qw(no_plan);

BEGIN { use_ok( 'Ends::NeighborMap' ); }
require_ok( 'Ends::NeighborMap' );
BEGIN { use_ok( 'Tree::Range' ); }
require_ok( 'Tree::Range' );

my $file = 't/test1.gff';

sub setup_files : Test(setup) {
}

sub test1 : Tests {
    ok(1,'meow');

    my ($numbins, $id_list_href, $trees_href) 
    = Ends::NeighborMap::make_trees(
        file => $file,
        prime => 5,
        flag => 0,
        distance => 5000,
        flag_6_distance => 5000,
        locus_tag => 'ID',
        binwidth => 100,
    );
    #say Dumper $trees_href;
    #say Dumper $trees_href->{'CHR1'}->search(8750,8751);
    #Chr1	TAIR8	gene	6790	8737	.	-	.	ID=AT1G01020;Name=AT1G01020;Note=ARV1
}

ends::make::trees->runtests;
