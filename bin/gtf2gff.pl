#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
END {close STDOUT}
$| = 1;

use FindBin;
use lib "$FindBin::Bin/../lib";
use GFF::Parser;

my $p = GFF::Parser->new(file => \*ARGV);

while (defined(my $gtf = $p->next())){
    $gtf->attribute_string(
        join ";", 
        map {
            my ($k, $v) = m/([^\s]+)\s"([^"]+)"/;
            defined $k ? (qq{$k=$v}) : ();
        } split /;\s/, $gtf->attribute_string()
    );
    say $gtf;
}
