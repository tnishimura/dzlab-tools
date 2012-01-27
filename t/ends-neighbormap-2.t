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
    my ($nm, $id, $strand, $overlapped, $focus,$minus, $plus) = @_;
    is_deeply([($nm->neighborhood($id))[0..4]],  [$focus, $strand, $overlapped, $minus, $plus], $id); 
}

#######################################################################
# Flag 0

sub flag_0_prime_5 : Tests {
    my $nm = Ends::NeighborMap->new(flag => 0, prime => 5);
    $nm->add('a','-', 9000, 19000);
    $nm->add('b','+', 10000,20000);
    $nm->add('c','-', 11000,21000);
    $nm->finalize();
    confirm($nm, 'b', '+', 1, 10000, 5000,  15000);
    confirm($nm, 'a', '-', 1, 19000, 14000, 24000);
}

sub flag_0_prime_3 : Tests {
    my $nm = Ends::NeighborMap->new(flag => 0, prime => 3);
    $nm->add('a','-', 9000, 19000);
    $nm->add('b','+', 10000,20000);
    $nm->add('c','-', 11000,21000);
    $nm->finalize();
    confirm($nm, 'b', '+', 1, 20000, 15000,  25000);
    confirm($nm, 'a', '-', 0, 9000, 4000, 14000);
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

    confirm($nm, 'target1', '+', 0, 10000, 8765, 12345);
    confirm($nm, 'target2', '-', 0, 21000, 19999, 26000);
    confirm($nm, 'target3', '+', 0, 1000, 888, 1050);
    confirm($nm, 'lonely', '-', 0, 110000, 105000, 115000);
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

    confirm($nm, 'target1', '+', 1, 14001, 10000, 19000);
    confirm($nm, 'target2', '-', 1, 19000, 14001, 21000);
    confirm($nm, 'target3', '+', 0, 1050, 1000, 3000);
    confirm($nm, 'lonely', '-', 0, 100000, 95000, 105000);
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

    confirm($nm, 'target1', '+', 1, 10000, 5000, 10350);
    confirm($nm, 'target2', '+', 0, 9000, 4000, 9999);
    confirm($nm, 'target3', '-', 0, 22000, 20002, 27000);
    confirm($nm, 'lonely', '-', 0, 110000, 105000, 115000);
}

sub flag_6_prime_3 : Tests {
    my $nm = Ends::NeighborMap->new(flag => 6, prime => 3);
    $nm->add('target1', '-', 10000, 10500);
    $nm->add('d1', '+', 10400, 20000);
    $nm->add('target2', '-', 9000, 10149);
    $nm->add('lonely','-', 100000, 110000);

    $nm->finalize();

    confirm($nm, 'target1', '-', 1, 10000, 5000, 10350);
    confirm($nm, 'target2', '-', 0, 9000, 4000, 9999);
    confirm($nm, 'lonely', '-', 0, 100000, 95000, 105000);
}

#######################################################################
# Flag 6

sub flag_7_prime_5 : Tests {
    my $nm = Ends::NeighborMap->new(flag => 7, prime => 5);
    #$nm->add('target1', '+', 10000, 10500);
    #$nm->add('d1', '+', 10400, 20000);
    #$nm->add('target2', '+', 9000, 10149);
    #$nm->add('target3', '-', 19852, 22000);
    $nm->add('lonely1','+', 100000, 110000);
    $nm->add('lonely2','-', 200000, 210000);

    $nm->finalize();
    #say Dumper $nm;

    #confirm($nm, 'target1', '+', 1, 10000, 5000, 10350);
    #confirm($nm, 'target2', '+', 0, 9000, 4000, 9999);
    #confirm($nm, 'target3', '-', 0, 22000, 20002, 27000);
    confirm($nm, 'lonely1', '+', 0, 100000, 90000, 110000);
    confirm($nm, 'lonely2', '-', 0, 210000, 200000, 220000);

    #my ($focus, $strand, $overlapped, $minus, $plus) = 
    #$nm->neighborhood('lonely2');
    #say Dumper {
    #    focus      => $focus ,
    #    strand     => $strand ,
    #    overlapped => $overlapped ,
    #    minus      => $minus ,
    #    plus       => $plus,
    #};
}

ends::neighborhood::Test->runtests;
