package binarysearch::test;
use base qw(Test::Class);
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More tests => 315;

sub setup_files : Test(setup) {
    my $self = shift;
    $self->{simple} = [ 17, 27, 28, 29, 35, 44, 49, 64, 65, 70, 83, 85, 89, 96, 98 ];
    $self->{randgen} = [ sort { $a<=> $b} map { -5000 + int rand 10000 } (0 .. 1000)];
}

sub extended_cmp : Tests{
    is(BinarySearch::_extended_cmp(3,5), -1);
    is(BinarySearch::_extended_cmp(4,3), 1);
    is(BinarySearch::_extended_cmp(3,3), 0);
    is(BinarySearch::_extended_cmp('inf',3), 1);
    is(BinarySearch::_extended_cmp(3,'-inf'), 1);
    is(BinarySearch::_extended_cmp(3,'inf'), -1);
    is(BinarySearch::_extended_cmp('inf', '-inf'), 1, '+inf vs -inf');
    is(BinarySearch::_extended_cmp('-inf', 'inf'), -1, '-inf vs +inf');
    is(BinarySearch::_extended_cmp('-inf', '-inf'), undef, '-inf vs -inf');
    is(BinarySearch::_extended_cmp('+inf', '+inf'), undef, '+inf vs +inf');
    is(BinarySearch::_extended_ge(4,3), 1);
    is(BinarySearch::_extended_ge(3,3), 1);
    is(BinarySearch::_extended_gt(4,3), 1);

    is(BinarySearch::_extended_le(3,3), 1);
    is(BinarySearch::_extended_le(3,4), 1);
    is(BinarySearch::_extended_lt(3,4), 1);
}

sub sanity : Tests {
    my @simple = @{shift->{simple}};
    is(15, scalar(@simple), "length");
}

sub linear : Tests {
    my @simple = @{shift->{simple}};
    for (@simple){
        is(BinarySearch::_greatest_lower_linear(\@simple, $_+1), $_, "_greatest_lower_linear");
        is(BinarySearch::_least_upper_linear(\@simple, $_-1), $_, "linear_upper_linear");
    }
    is(BinarySearch::_greatest_lower_linear(\@simple, 0), undef, "_greatest_lower_linear");
    is(BinarySearch::_least_upper_linear(\@simple, 100), undef, "_least_upper_linear");
}

sub binary : Tests {
    my @simple = @{shift->{simple}};
    for (@simple){
        is(greatest_lower(\@simple, $_+1), $_, "greatest_lower");
        is(least_upper(\@simple, $_-1), $_, "least_upper");
    }
    is(greatest_lower(\@simple, 0), undef, "great_lower inf");
    is(least_upper(\@simple, 100), undef, "least_upper inf");

    for (@simple){
        is(
            BinarySearch::_greatest_lower_linear(\@simple, $_+1), 
            greatest_lower(\@simple, $_+1), 
            "linear vs binary"
        );
        is(
            BinarySearch::_least_upper_linear(\@simple, $_-1), 
            least_upper(\@simple, $_-1), 
            "linear vs binary"
        );
    }
}
sub binary_rand : Tests {
    my @randarr = @{shift->{randgen}};

    for (map { -6000 + int rand 12000 } (0 .. 100)){
        is(
            BinarySearch::_greatest_lower_linear(\@randarr, $_+1), 
            greatest_lower(\@randarr, $_+1), 
            "linear vs binary gl"
        );
        is(
            BinarySearch::_least_upper_linear(\@randarr, $_-1), 
            least_upper(\@randarr, $_-1), 
            "linear vs binary lu"
        );
    }
}


BEGIN { 
    use_ok( 
        'BinarySearch', 
        qw/greatest_lower least_upper/
    ); 
}
require_ok( 'BinarySearch' );

binarysearch::test->runtests;
