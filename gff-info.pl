#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Statistics qw/gff_info/;
use YAML qw/Dump/;

END {close STDOUT}
$| = 1;

for my $file (@ARGV) {
    print Dump gff_info $file;
}

