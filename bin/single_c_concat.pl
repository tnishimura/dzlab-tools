#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib";
use autodie;
use Config::General qw(ParseConfig);
use Data::Dumper;
use feature 'say';
use DZUtil qw/single_c_concat/;

die "usage: $0 contextual-single-c-files.gff ..." unless @ARGV;

single_c_concat(@ARGV);
