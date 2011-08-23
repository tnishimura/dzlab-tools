#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;

END {close STDOUT}
$| = 1;

if (@ARGV != 2){
    say "usage: $0 superset supposed-subset";
    exit 1;
}

my $superset_parser = GFF::Parser->new(file => $ARGV[0]);
my $subset_parser = GFF::Parser->new(file => $ARGV[1]);

my %superset_positions;

while (defined(my $gff = $superset_parser->next())){
    $superset_positions{$gff->sequence,$gff->start} = $gff->score;
}

while (defined(my $gff = $subset_parser->next())){
    if ($gff->score > 0 && ! exists $superset_positions{$gff->sequence,$gff->start}){
        say $gff->to_string();
    }
}
