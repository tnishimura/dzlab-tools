#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;

END {close STDOUT}
$| = 1;

while (defined(my $line = <ARGV>)){
    chomp $line;
    my ($bsread, undef, $read) = map {scalar <ARGV>} (1..3);
    compare($read,$bsread);
}
print_total();

{
    my $total = 0;
    sub compare{
        my ($read_str, $bsread_str) = @_;
        my @read = split //, $read_str;
        my @bsread = split //, $bsread_str;
        if (@read != @bsread){
            die "$. read len don't match";
        }
        for (0 .. $#read){
            my ($r, $b) = ($read[$_], $bsread[$_]);
            if ($r ne $b){
                if ($r eq 'C' && $b eq 'T'){
                    $total++;

                }
                else {
                    die "non C->T change";
                }
            }
        }
    }
    sub print_total{
        say "total: $total";
    }
}
