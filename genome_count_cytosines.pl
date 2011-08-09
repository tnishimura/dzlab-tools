#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;

END {close STDOUT}
$| = 1;

my $count = 0;
while (defined(my $line = <ARGV>)){
    chomp $line;
    next if $line =~ /^>/;
    my @cg = $line =~ /C|G/g;
    $count += scalar @cg;
}
say $count;
