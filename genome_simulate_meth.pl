#!/usr/bin/env perl
use 5.10.0;
use strict;
use warnings FATAL => "all";
use Data::Dumper;
use feature 'say';
use autodie;
use Pod::Usage;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/lib";
use FastaReader;

END {close STDOUT}

my $bs_rate = 0.25;
my $prefix;

my $result = GetOptions (
    "bs-rate|b=f" => \$bs_rate,
    "prefix|p=s"  => \$prefix,
);

unless (@ARGV == 1 && -f $ARGV[0]){
    say "usage: $0 [-b bs-rate] [-p output_prefix] input";
    exit 1;
}
my $input = $ARGV[0];

if (! defined $prefix){ 
    $prefix = $input;
}

open my $coord_fh, '>', ($prefix . ".simmeth.gff");
open my $output_fh, '>', ($prefix . ".simmeth");

my $fr = FastaReader->new(file => $input, slurp => 1);

my %iterators = $fr->make_iterators();
say Dumper \%iterators;

while (my ($seq,$iter) = each %iterators) {
    say $output_fh ">$seq";
    my $position = 1;
    
    while (defined(my $line = $iter->())){
        my @bases = split //, $line;
        BASE:
        for my $index (0 .. $#bases){
            my $base = uc $bases[$index];
            if ($base eq 'C'){
                my $context = $fr->get_context($seq, $position, base => 1);

                if (rand(1) < $bs_rate){
                    say $coord_fh join "\t", $seq, "C", $context, $position, $position, qw/. + . ./;
                }
                else {
                    say $coord_fh join "\t", $seq, "T", $context, $position, $position, qw/. + . ./;
                    $bases[$index] = 'T';
                }
            }

            ++$position;
        }

        say $output_fh join "", @bases;
    }
}
close $output_fh;
close $coord_fh;

