#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;

use FindBin;
use lib "$FindBin::Bin/lib";
use DZUtil qw/fastq_read_length/;
use Getopt::Long;

my $result = GetOptions (
    "no-filename|n" => \(my $no_filename),
);

if (!$result || @ARGV == 0){
    say "usage: $0 alignment.fastq .. ";
    exit 1;
}

for my $file (@ARGV) {
    # body...
    my $length = fastq_read_length($file);

    if (! defined $length){
        say "$0: fastq file $file malformed?";
    }
    if ($no_filename){
        say "$length";
    }
    else{
        say "$file: $length";
    }
}

