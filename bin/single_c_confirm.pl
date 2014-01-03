#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/../lib";
use FastaReader;
use GFF::Parser;
use Pod::Usage;
use Getopt::Long;

END {close STDOUT}
$| = 1;

my $result = GetOptions (
    "reference|r=s" => \(my $reference),
);

if (!$result || ! $reference || ! -f $reference || grep { ! -f } @ARGV){
    say "check that single-c file positions are actually C's.";
    say "usage: single_c_confirm.pl -r reference single-c.gff ...";
    exit 1;
}

my $fr = FastaReader->new(file => $reference, slurp => 1);

for my $file (@ARGV) {
    my $p = GFF::Parser->new(file => $file, normalize => 0);
    while (defined(my $gff = $p->next())){
        my $c = $gff->get_column('c');
        my $seq = $gff->sequence;
        my $pos = $gff->start;
        my $strand = $gff->strand;

        if (defined $c && $c > 0){
            my $base = $fr->get($seq, $pos, $pos);
            if ($base ne 'C' && $base ne 'G'){
                die("error at $seq $pos: not C/G");
            }
            if (defined $strand){
                if ($strand eq '+' && $base ne 'C'){
                    die("error at $seq $pos + ");
                }
                elsif ($strand eq '-' && $base ne 'G'){
                    die("error at $seq $pos -");
                }
            }
        }
    }
}
