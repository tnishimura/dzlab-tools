#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use List::MoreUtils qw/all/;
use Term::ProgressBar;

for my $file (@ARGV) {

    my $size = (stat($file))[7];

    say "$file:";

    my $pb = Term::ProgressBar->new({count => $size});
    my $counter = 0;
    $pb->minor(0);

    open my $in, '<', $file;

    while (! eof $in){
        my @lines = map { scalar <$in> } (1 .. 4);
        my @lens = map { length $_ } @lines;

        # all undef
        if (all { ! defined $_ } @lines){
            last;
        }
        else{
            # only some undef
            if (grep { ! defined $_ } @lines){
                die "uneven number of lines @ around $."
            }
            elsif ($lines[0] !~ /^@/ 
                || $lines[2] !~ /^\+/ 
                # third line could be just a lone '+', in newer fastq files
                || ($lens[0] != $lens[2] && $lines[2] !~ /^\+$/)
                || $lens[1] != $lens[3] 
                ){
                die "malformed FASTQ quad @ $.\n" . join "", @lines;
            }
        }
        $pb->update(tell($in)) if ++$counter % 100_000 == 0;
    }
    $pb->update($size);

    close $in;
}
say "OK";
exit 0;
