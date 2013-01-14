#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Scalar::Util qw/looks_like_number/;

use Pod::Usage;
use Getopt::Long;

if (! @ARGV && -t STDIN){
    pod2usage(-verbose => 2, -noperldoc => 1);
}

while (defined(my $line = <ARGV>)){
    chomp $line;

    my ($qname, $flag, $rname, $leftmost, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual) 
    = split /\t/, $line;

    if ($line =~ /^@/ or (looks_like_number $flag && ! ($flag & 0x4))) {
        say $line;
    }
}

=head1 NAME

sam-only-mapped.pl in.sam > out.sam

=cut

