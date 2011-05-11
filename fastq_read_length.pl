#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;

use FindBin;
use lib "$FindBin::Bin/lib";
use DZUtil qw/fastq_read_length/;

if (@ARGV != 1){
    say "usage: $0 alignment.fastq";
    exit 1;
}

my $length = fastq_read_length($ARGV[0]);

if (! defined $length){
    say "$0: fastq file malformed?";
}
say $length;
