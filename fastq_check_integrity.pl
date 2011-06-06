#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use List::MoreUtils qw/all/;

my $counter = 0;
while (1){
    my @lines = map { scalar <> } (1 .. 4);
    # all undef
    if (all { ! defined $_ } @lines){
        last;
    }
    else{
        # only some undef
        if (grep { ! defined $_ } @lines){
            die "uneven number of lines @ around $."
        }
        elsif ($lines[0] !~ /^@/ || $lines[2] !~ /^\+/ || length $lines[0] != length $lines[2] ||  length $lines[1] != length $lines[3]){
            die "malformed FASTQ quad @ $.\n" . join "", @lines;
        }
    }

    say $counter if --$counter % 100_000 == 0;
}
say "OK";
