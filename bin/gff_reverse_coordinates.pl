#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/../lib";
use FastaReader;


END {close STDOUT}
$| = 1;

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "reference|r=s" => \(my $reference),
    "input|i=s" => \(my $input),
);
if (!$results || ! defined $input || ! defined $reference){
    say "usage: $0 --reference original-genome.fasta --input input.gff";
    exit 1;
}


{
    my $fr = FastaReader->new(file => $reference);
    my %seqlen = $fr->sequence_lengths();

    sub coord_reverse{
        my ($seq, $start, $end) = @_;
        my $len = $seqlen{$seq};
        return ($len - $end + 1, $len - $start + 1);
    }
}
sub reverse_strand{
    given (shift){
        when (undef) { return '.' }
        when ('+') { return '-' }
        when ('-') { return '+' }
        default { die "huh?"; }
    }
}

open my $gff_fh, '<', $input;
while (defined(my $gff_line = <$gff_fh>)){
    chomp $gff_line;
    my @gff = split /\t/, $gff_line;
    if (@gff != 9){
        say $gff_line;
        next;
    }
    @gff[3,4] = coord_reverse(@gff[0,3,4]);
    $gff[6] = reverse_strand($gff[6]);
    say join "\t", @gff;
}
close $gff_fh;
