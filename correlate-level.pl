#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;

END {close STDOUT}
$| = 1;

use Pod::Usage;
use Getopt::Long;

my $help;
my $verbose;
my $result = GetOptions (
    "verbose"  => \$verbose,
    "help"     => \$help,
    "ct=s"     => \my $ct_cutoff_arg,
    "methyl=s" => \my $methyl_cutoff_arg,
);

my $file = shift;
if (!$result || $help || ! $ct_cutoff_arg || ! $methyl_cutoff_arg || ! $file){
    say "$0: take correlate*Ends.pl output and stratify by ct and methyl upper limit cutoffs";
    say "usage: $0 --ct 2,4,5,7 --methyl 10,20,50,80 file.gff";
    exit 1;
}

my @ct_cutoffs = split /,/, $ct_cutoff_arg;
my @methyl_cutoffs = split /,/, $methyl_cutoff_arg;

my %fh;
for my $ct (@ct_cutoffs) {
    for my $methyl (@methyl_cutoffs) {
        open my $output, '>', sprintf("%s.%d_%03d", $file, $ct, $methyl);
        $fh{$ct}{$methyl} = $output;
    }
}

open my $input, '<', $file;
while (defined(my $line = <$input>)){
    chomp $line;
    my @split = $line;
    my ($chh, $thh) = @split[13,14];
    my $ct = $chh + $thh;
    my $original = join "\t", @split[0..8];

    for my $ct_cutoff (@ct_cutoffs) {
        for my $methyl_cutoff (@methyl_cutoffs) {
            my $output = $fh{$ct_cutoff}{$methyl_cutoff};
            if ($ct == 0){
                say $output $original;
            }
            else {
                my $methyl = $chh/$ct;
                if ($methyl*100 <= $methyl_cutoff && $ct <= $ct_cutoff){
                    say $output $original;
                }
            }
        }
    }
}
close $input;
close $_ for (values %fh);
