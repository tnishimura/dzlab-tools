#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;

if (@ARGV != 2){
    say "usage: $0 <seqa> <seqb>";
    exit 255;
}

if (length $ARGV[0] != length $ARGV[1]){
    say "<seqa> and <seqb> should be same length";
    exit 255;
}
my ($seqa, $seqb) = @ARGV;
$seqa = uc $seqa;
$seqb = uc $seqb;

my @splita = split //, $seqa;
my @splitb = split //, $seqb;

print "\n";
say $seqa;
while ((my $aa = shift @splita) && (my $bb = shift @splitb)){
    if ($aa eq $bb){
        print "-";
    } else {
        print "|";
    }
}
print "\n";
say $seqb;

for my $i (1 .. length $seqa) {
    if ($i % 5 == 0){
        print "    |";
    }
}
print "\n";

for my $i (1 .. length $seqa) {
    if ($i % 5 == 0){
        printf "%5d", $i;
    }
}
print "\n";


=head1 NAME

seqdiff.pl - Visually compare two sequences of identical length.

=head1 SYNOPSIS

Usage examples:

 seqdiff.pl acccct aCCCCt

Output:
 
 ACCCCT
 ------
 ACCCCT
     |
     5

=cut

