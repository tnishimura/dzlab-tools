#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
END {close STDOUT}
$| = 1;
use FindBin;
use lib "$FindBin::Bin/lib";
use CHI::Memo qw/dump clear_cache/;

if (@ARGV != 1){
    say "$0 clear";
    say "$0 dump";
}

for ($ARGV[0]) {
    when ('clear'){
        clear_cache();
    }
    when ('dump'){
        dump();
    }
}


