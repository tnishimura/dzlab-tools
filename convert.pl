#!/usr/bin/env perl
use strict;
use warnings;
use Pod::Usage;

if (@ARGV != 2 || ($ARGV[0] ne 'c2t' && $ARGV[0] ne 'g2a')){
    pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]);
}

my $pattern   = $ARGV[0];
my $fastafile = $ARGV[1];

open(FAS,"<$fastafile") or die "Can't read input file";
while(my $line = <FAS>) {
    if($line =~ m/^[ACGTN]+/i) {
        $line =~ tr/Cc/Tt/ if $pattern eq 'c2t';
        $line =~ tr/Gg/Aa/ if $pattern eq 'g2a';
    }
    print $line;
}
close(FAS);


=head1 NAME

convert.pl - Convert a fasta file (either c->t or g->a).

=head1 SYNOPSIS

Usage examples:

 convert.pl c2t genome.fasta > genome.c2t.fasta
 convert.pl g2a genome.fasta > genome.g2a.fasta

=back

=cut
