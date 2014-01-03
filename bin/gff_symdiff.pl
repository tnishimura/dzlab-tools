#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/../lib";
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

=head1 NAME

gff_symdiff.pl - given a list of GFF files (say a.gff, b.gff, c.gff), will
create files with the .UNIQ prefix (a.gff.UNIQ, b.gff.UNIQ, c.gff.UNIQ) which
contain entries that are unique to each.  Equality is determined by comparing
the sequence, start and end coordinate.  

=head1 SYNOPSIS

Usage examples:

 gff_symdiff.pl a.gff b.gff c.gff

=cut

