package ends::neighborhood::Test;
use base qw(Test::Class);
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More qw(no_plan);
use Ends::NeighborMap;
use List::Util qw/min max/;

sub confirm{
    my ($nm, $id, $strand, $focus,$minus, $plus) = @_;
    is_deeply([$nm->neighborhood($id)],  [$focus, $strand, $minus, $plus]); 
}

#######################################################################
# Flag 0

sub flag_0_prime_5 : Tests {
    my $nm = Ends::NeighborMap->new(flag => 0, prime => 5);
    $nm->add('a','-', 9000, 19000);
    $nm->add('b','+', 10000,20000);
    $nm->add('c','-', 11000,21000);
    $nm->finalize();
    confirm($nm, 'b', '+', 10000, 5000,  15000);
    confirm($nm, 'a', '-', 19000, 14000, 24000);
}

sub flag_0_prime_3 : Tests {
    my $nm = Ends::NeighborMap->new(flag => 0, prime => 3);
    $nm->add('a','-', 9000, 19000);
    $nm->add('b','+', 10000,20000);
    $nm->add('c','-', 11000,21000);
    $nm->finalize();
    confirm($nm, 'b', '+', 20000, 15000,  25000);
    confirm($nm, 'a', '-', 9000, 4000, 14000);
}

#######################################################################
# Flag 2

sub flag_2_prime_5 : Tests {
    my $nm = Ends::NeighborMap->new(flag => 2, prime => 5);
    $nm->add('d1','-', 12345, 19999);
    $nm->add('u1','-', 3000, 8765);
    $nm->add('u2','-', 500, 888);
    $nm->add('target1','+', 10000, 14001);
    $nm->add('target2','-', 19000, 21000);
    $nm->add('target3','+', 1000, 1050);
    $nm->add('lonely','-', 100000, 110000);
    $nm->finalize();

    confirm($nm, 'target1', '+', 10000, 8765, 12345);
    confirm($nm, 'target2', '-', 21000, 19999, 26000);
    confirm($nm, 'target3', '+', 1000, 888, 1050);
    confirm($nm, 'lonely', '-', 110000, 105000, 115000);
}

sub flag_2_prime_3 : Tests {
    my $nm = Ends::NeighborMap->new(flag => 2, prime => 3);
    $nm->add('d1','-', 12345, 19999);
    $nm->add('u1','-', 3000, 8765);
    $nm->add('u2','-', 500, 888);
    $nm->add('target1','+', 10000, 14001);
    $nm->add('target2','-', 19000, 21000);
    $nm->add('target3','+', 1000, 1050);
    $nm->add('lonely','-', 100000, 110000);
    $nm->finalize();

    confirm($nm, 'target1', '+', 14001, 10000, 19000);
    confirm($nm, 'target2', '-', 19000, 14001, 21000);
    confirm($nm, 'target3', '+', 1050, 1000, 3000);
    confirm($nm, 'lonely', '-', 100000, 95000, 105000);
}

#######################################################################
# Flag 6

sub flag_6_prime_5 : Tests {
    my $nm = Ends::NeighborMap->new(flag => 6, prime => 5);
    $nm->add('target1', '+', 10000, 10500);
    $nm->add('d1', '+', 10400, 20000);
    $nm->add('target2', '+', 9000, 10149);
    $nm->add('target3', '-', 19852, 22000);
    $nm->add('lonely','-', 100000, 110000);

    $nm->finalize();

    confirm($nm, 'target1', '+', 10000, 5000, 10350);
    confirm($nm, 'target2', '+', 9000, 4000, 9999);
    confirm($nm, 'target3', '-', 22000, 20002, 27000);
    confirm($nm, 'lonely', '-', 110000, 105000, 115000);
}

sub flag_6_prime_3 : Tests {
    my $nm = Ends::NeighborMap->new(flag => 6, prime => 3);
    $nm->add('target1', '-', 10000, 10500);
    $nm->add('d1', '+', 10400, 20000);
    $nm->add('target2', '-', 9000, 10149);
    $nm->add('lonely','-', 100000, 110000);

    $nm->finalize();

    confirm($nm, 'target1', '-', 10000, 5000, 10350);
    confirm($nm, 'target2', '-', 9000, 4000, 9999);
    confirm($nm, 'lonely', '-', 100000, 95000, 105000);
}

ends::neighborhood::Test->runtests;
