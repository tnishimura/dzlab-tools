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
use Launch;

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    #"longname|l=isf" => \my $var,
    "reference|r=s" => \my $reference,
    "bowtie|b=s" => \my $bowtie,
    "window-by-fixed|w=i" => \(my $window_size = 50),
);
if (!$result || !$bowtie || !$reference){
    say "usage: bowtie_window.pl -b bowtie -r genome";
    exit 1;
}

my $basename = $bowtie;
$basename =~ s/\.bowtie$//;

my $gff = "$basename.gff";
my $w50 = "$basename.w$window_size.gff";

launch("perl -S parse_bowtie.pl -g --windowable -o ?? $bowtie",
    expected => $gff,
);
launch("perl -S window_by_fixed.pl --count-in-scores -n -w $window_size -k -r $reference -o ?? $gff",
    expected => $w50,
);
