#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use List::Util qw/sum/;
use Statistics::Descriptive;

END {close STDOUT}
$| = 1;

use YAML qw/Bless Dump/;

my %yaml;

for my $file (@ARGV) {
    my $count = 0;
    my $lengths = Statistics::Descriptive::Full->new();
    my %file_yaml;

    my $parser = GFF::Parser->new(file => $file, normalize => 0);
    while (defined(my $gff = $parser->next())){
        my ($seq, $feature) = ($gff->sequence(), $gff->feature());

        $file_yaml{sequences}{$seq}++;
        $file_yaml{features}{$feature}++;
        $file_yaml{count}++;

        $lengths->add_data($gff->end - $gff->start + 1);
    }

    $file_yaml{mean}                = sprintf "%.1f", $lengths->mean();
    $file_yaml{standard_deviation}  = sprintf "%.1f", $lengths->standard_deviation();
    $file_yaml{min}                 = $lengths->min();
    $file_yaml{median}              = $lengths->median();
    $file_yaml{max}                 = $lengths->max();
    $file_yaml{number_of_sequences} = scalar(keys %{$file_yaml{sequences}});
    $file_yaml{number_of_features}  = scalar(keys %{$file_yaml{features}});

    Bless(\%file_yaml)->keys([qw/
        count mean standard_deviation min median max number_of_sequences sequences number_of_features features
        /]);

    $yaml{$file} = \%file_yaml;
}

say Dump \%yaml;
