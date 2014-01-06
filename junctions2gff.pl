#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;

use Pod::Usage;
use Getopt::Long;

my $label;
my $result = GetOptions (
    "label|l=s" => \$label,
);
if (!$label){
    say "usage: $0 -l s_4 < junctions.bed > output.gff";
    exit 1;
}

<>; #skip header

while (defined(my $line = <>)){
    chomp $line;
    my ($seq, $start, $end, $juncname, $score, $strand, undef, undef, undef, undef, $blocksizes, $blockstarts) 
    = split /\t/, $line;

    say join "\t", $seq, $juncname, "$label.junction", $start, $end, $score, $strand, '.', "ID=$juncname";

    my @sizes = split /,/, $blocksizes;
    my @starts = split /,/, $blockstarts;

    if (@sizes != 2 || @starts != 2){
        die "$line: block count != 2?";
    }

    my $counter = 0;

    for my $i (0, 1) {
        say join "\t", $seq, $juncname, "$label.block", $start + $starts[$i], $start + $starts[$i]+$sizes[$i]-1, $score, $strand, '.', "ID=$juncname.$i";
        # body...
    }
}

__DATA__
==> s_5.gff <==
chrc    JUNC00000001    s5.junction     97695   106949  1       -       .       ID=JUNC00000001
chrc    JUNC00000001    s5.block        97695   97703   1       -       .       ID=JUNC00000001.0
chrc    JUNC00000001    s5.block        106922  106948  1       -       .       ID=JUNC00000001.1
chrc    JUNC00000002    s5.junction     131699  140953  1       +       .       ID=JUNC00000002
chrc    JUNC00000002    s5.block        131699  131725  1       +       .       ID=JUNC00000002.0
chrc    JUNC00000002    s5.block        140944  140952  1       +       .       ID=JUNC00000002.1
chr2    JUNC00000003    s5.junction     17344   23720   21      -       .       ID=JUNC00000003
chr2    JUNC00000003    s5.block        17344   17361   21      -       .       ID=JUNC00000003.0
chr2    JUNC00000003    s5.block        23692   23719   21      -       .       ID=JUNC00000003.1
chr2    JUNC00000004    s5.junction     17344   29110   21      -       .       ID=JUNC00000004

==> junctions.bed <==
track name=junctions description="TopHat junctions"
chrc    97695   106949  JUNC00000001    1       -       97695   106949  255,0,0 2       9,27    0,9227
chrc    131699  140953  JUNC00000002    1       +       131699  140953  255,0,0 2       27,9    0,9245
chr2    17344   23720   JUNC00000003    21      -       17344   23720   255,0,0 2       18,28   0,6348
chr2    17344   29110   JUNC00000004    21      -       17344   29110   255,0,0 2       18,28   0,11738
chr2    23704   29110   JUNC00000005    21      -       23704   29110   255,0,0 2       18,28   0,5378
chr2    55908   74536   JUNC00000006    7       -       55908   74536   255,0,0 2       8,28    0,18600
chr2    74553   85675   JUNC00000007    4       -       74553   85675   255,0,0 2       29,10   0,11112
chr2    74553   74816   JUNC00000008    39      -       74553   74816   255,0,0 2       32,32   0,231
chr2    75751   75878   JUNC00000009    4       +       75751   75878   255,0,0 2       29,21   0,106

