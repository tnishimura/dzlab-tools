#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Pod::Usage;
use Getopt::Long;

if (! @ARGV && -t STDIN){
    pod2usage(-verbose => 2, -noperldoc => 1);
}

while (<ARGV>){
    my ($readid, undef, undef, undef, $seq, $qual) = split /\t/, $_;
    print("\@$readid\n$seq\n+$readid\n$qual\n");
}

=head1 bowtie2fastq.pl 

Recover FASTQ from Bowtie...

 bowtie2fastq.pl < in.bowtie > out.fastq

=cut
