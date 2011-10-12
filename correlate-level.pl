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
    "chh=s"     => \my $chh_cutoff_arg,
    "methyl=s" => \my $methyl_cutoff_arg,
);

my $file = shift;
if (!$result || $help || ! $chh_cutoff_arg || ! $methyl_cutoff_arg || ! $file){
    say "$0: take correlate*Ends.pl output and stratify by chh and methyl upper limit cutoffs";
    say "usage: $0 --chh 2,4,5,7 --methyl 10,20,50,80 file.gff";
    exit 1;
}

my @chh_cutoffs = split /,/, $chh_cutoff_arg;
my @methyl_cutoffs = split /,/, $methyl_cutoff_arg;

my %fh;
for my $chh (@chh_cutoffs) {
    for my $methyl (@methyl_cutoffs) {
        open my $output, '>', sprintf("%s.%d_%03d", $file, $chh, $methyl);
        $fh{$chh}{$methyl} = $output;
    }
}

open my $input, '<', $file;
while (defined(my $line = <$input>)){
    say STDERR $. if $. % 50000 == 0;
    chomp $line;
    my @split = split /\t/, $line;
    my ($chh, $thh) = @split[13,14];
    my $ct = $chh + $thh;
    my $original = join "\t", @split[0..8];

    for my $chh_cutoff (@chh_cutoffs) {
        for my $methyl_cutoff (@methyl_cutoffs) {
            my $output = $fh{$chh_cutoff}{$methyl_cutoff};
            if ($ct == 0){
                say $output $original;
            }
            else {
                my $methyl = $chh/$ct;
                if ($methyl*100 <= $methyl_cutoff && $chh <= $chh_cutoff){
                    say $output $original;
                }
            }
        }
    }
}
close $input;
