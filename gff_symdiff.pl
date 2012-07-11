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

my @outfh;
my %set;

while (my ($index,$file) = each @ARGV) {
    my $p = GFF::Parser->new(file => $file, normalize => 1);
    while (defined(my $gff = $p->next())){
        push @{$set{$gff->sequence // '.'}{$gff->start}{$gff->end}}, [$index, $gff->to_string];
    }
    my $outfile = $file;
    open my $fh, '>', $file . ".UNIQ";

    push @outfh, $fh;
}

for my $seqid (sort keys %set) {
    for my $start (sort keys %{$set{$seqid}}) {
        for my $end (sort keys %{$set{$seqid}{$start}}) {
            if (@{$set{$seqid}{$start}{$end}} == 1){
                my $igff = $set{$seqid}{$start}{$end}->[0];
                say {$outfh[$igff->[0]]} $igff->[1];
            }
        }
    }
}
