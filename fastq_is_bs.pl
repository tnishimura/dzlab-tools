#!/usr/bin/env perl
use 5.10.0;
use strict;
use warnings FATAL => "all";
use Data::Dumper;
use feature 'say';
use autodie;

END {close STDOUT}

for my $file (@ARGV) {
    open my $fh, '<', $file;
    my ($a, $c, $g, $t) = (0) x 4;

    my $length;

    READ:
    while (defined(my $line = <$fh>)){
        chomp $line;
        last READ unless 1 .. 100_000;
        next READ unless $. % 4 == 2;
        $length //= length $line;
        if (length $line != $length){
            die "$file: corrupt?";
        }
        $line = lc $line;

        $a += $line =~ tr/a/a/;
        $c += $line =~ tr/c/c/;
        $g += $line =~ tr/g/g/;
        $t += $line =~ tr/t/t/;
    }

    my $total = $a + $c + $g + $t;
    my $c_ratio = $c / $total;

    printf "%s: %.4f %s\n", $file, $c_ratio, ($c_ratio < .125 ? 'BS' : 'gDNA');
    #printf "a: %.4f\n", $a / $total;
    #printf "c: %.4f\n", $c / $total;
    #printf "g: %.4f\n", $g / $total;
    #printf "t: %.4f\n", $t / $total;
}

