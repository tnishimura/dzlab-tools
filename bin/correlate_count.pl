#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;

END {close STDOUT}
$| = 1;

my ($nm,$match) = (0,0);
while (defined(my $line = <ARGV>)){
    chomp $line;
    if ($line =~ /^\./){
        $nm++;
    }
    else{
        $match++;
    }
}
printf("%d out of %d reads aligned\n", $match, $nm + $match);
