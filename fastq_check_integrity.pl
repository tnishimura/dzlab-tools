#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use List::MoreUtils qw/all/;

my $counter = 0;
for my $file (@ARGV) {

    my $size = (stat($file))[7];

    open my $in, '<', $file;

    while (! eof $in){
        my @lines = map { scalar <$in> } (1 .. 4);

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
        printf("%f\n",tell($in)/$size) if ++$counter % 100_000 == 0;
    }

    close $in;
}
say "OK";
