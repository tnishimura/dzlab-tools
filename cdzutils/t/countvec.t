#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use Data::Dumper;
use autodie;
use Test::More;
use Test::Exception;
use CDZUtils::CountVec;

my $cv = CDZUtils::CountVec->new(size => 100, base => 1, initial => 0);

$cv->increment_range(1, 100, 1);
is_deeply($cv->get_int(1, 100), [(1) x 100]);

for (2 .. 100){
    $cv->increment_range($_, 100, 1);
}

is_deeply($cv->get_int, [1 .. 100]);

done_testing();
