#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use FastaReader;


END {close STDOUT}
$| = 1;
use Pod::Usage;
use Getopt::Long;

FastaReader::rc_file($ARGV[0], \*STDOUT);
