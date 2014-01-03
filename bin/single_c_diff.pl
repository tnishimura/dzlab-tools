#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use YAML qw/Dump/;

END {close STDOUT}
$| = 1;

if (@ARGV < 2){
    say "compare genome_shear mock single-c to bs-sequel output";
    say "usage: $0 superset supposed-subset";
    exit 1;
}

my $superset_parser = GFF::Parser->new(file => shift);

my %superset_positions;
my %stats = (
    found_and_underscore => 0,
    correct => 0,
    overscore => 0,
    error => 0,
    total_real_c => 0,
    total_detected_c => 0,
    non_methylated_detected => 0,
);

while (defined(my $gff = $superset_parser->next())){
    $superset_positions{$gff->sequence,$gff->start} = $gff->score;
    $stats{total_real_c}+=$gff->score;
}

for my $singlec (@ARGV) {
    say STDERR "processing $singlec";
    my $subset_parser = GFF::Parser->new(file => $singlec);
    while (defined(my $gff = $subset_parser->next())){
        my $subset_count = $gff->score > 0 ? $gff->get_column('c') : 0;
        my $superset_count = 
        exists $superset_positions{$gff->sequence,$gff->start} ?  $superset_positions{$gff->sequence,$gff->start} : 0;

        if ($subset_count > 0 && $superset_count ==0){
            $stats{error}++;
        }
        elsif ($subset_count > 0 && $superset_count > $subset_count){
            $stats{found_and_underscore}++;
        }
        elsif ($subset_count > 0 && $superset_count == $subset_count){
            $stats{correct}++;
        }
        elsif ($subset_count > 0 && $superset_count < $subset_count){
            $stats{overscore}++;
        }
        elsif ($subset_count == 0){
            $stats{non_methylated_detected}++;
        }
        $stats{total_detected_c}+=$subset_count;
    }

}

say Dump \%stats;
