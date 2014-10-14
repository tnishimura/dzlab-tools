#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use v5.12.0;
use Data::Dumper;
use autodie;

use FindBin;
use lib "$FindBin::Bin/lib";
use BowtieParser;
use FastaReader;
use BigArray;

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "reference|r=s" => \(my $reference),
    "feature|f=s" => \(my $feature),
    # "max-entropy|m=i" => \(my $max_entropy = .17),
    "max-entropy|m=i" => \(my $max_entropy),
);
usage() if !$result || !$feature || ! $reference;

my $file = shift || usage();

my %stats;
my %position = ( A => 1, C => 2, G => 3, T => 4,);

my $fastareader = FastaReader->new(file => $reference, slurp => 0);
my %coverage = map { 
    uc($_) => BigArray->new(base => 1, size => $fastareader->get_length($_)) 
} $fastareader->sequence_list();

# general idea:
# record

my $p = BowtieParser->new(file => $file);
while (defined(my $alignment = $p->next(1))){
    my ($readid, $strand, $seqid, $start, $read, $quality, $mystery, $mismatches) = @$alignment;
    $coverage{uc($seqid)}->increment_range($start, $start + length($read) - 1);
    MM:
    for my $mm (@$mismatches) {
        my ($abs_coord, $base_in_ref, $base_in_read) = @$mm;
        next MM if ($base_in_read eq 'N');

        if (!$stats{$seqid}{$abs_coord}){
            $stats{$seqid}{$abs_coord} = [$base_in_ref, 0, 0, 0, 0];
        }
        my $record = $stats{$seqid}{$abs_coord};

        if ($record->[0] eq $base_in_ref){
            $record->[$position{$base_in_read}]++;
        }
        else{
            die "Bowtie error";
        }
        #say join "\t", $seqid, $abs_coord, $base_in_ref, $base_in_read;
    }
    counter();
}

use List::Util qw/max maxstr sum/;

for my $seq (keys %stats) {
    my $bigarray = $coverage{uc($seq)};

    POS:
    for my $pos (keys %{$stats{$seq}}) {
        my $cover = $bigarray->get($pos);
        my $record = $stats{$seq}{$pos}; # [$base_in_ref, $A, $C, $G, $T]
        my $base_in_ref = $record->[0];
        
        my $total_snps  = sum @{$record}[1 .. 4];
        my $stayed_same = $cover - $total_snps;

        die "bug" if ($record->[$position{$base_in_ref}] != 0);

        # ad-hoc: if the total number of snps is significantly less than non-snps, don't report
        next POS  if ($total_snps * 10 < $stayed_same); 

        # update record so the base_in_ref entry is number of reads reporting non-snps.
        $record->[$position{$base_in_ref}] = $stayed_same;

        my $entropy = entropy(@{$record}[1 .. 4]);
        if (! defined $max_entropy || $entropy < $max_entropy){
            # score is -log($entropy) so that higher scores = less entropy
            say join "\t", lc($seq), $base_in_ref, $feature, $pos, $pos, -log($entropy), qw/. ./, 
            "A=$record->[1];C=$record->[2];G=$record->[3];T=$record->[4];entropy=$entropy";
        }
    }
}

#######################################################################
# utility

sub entropy{
    my @nums = map { $_ + 1 } @_;
    my $total = sum @nums;
    my @freq  = map { $_ / $total } @nums;
    return (- (sum map { $_ == 0 ? 0 : $_ * log $_ } @freq));
}
sub usage{
    print STDERR <<'END';
 bowtie-extract-snps.pl - given a bowtie alignment file, create a list of detected SNPs over a certain entry

 bowtie-extract-snps.pl -f feature_name [-m .17] bowtie-alignment > snps.gff
END
    exit 1;
}

sub counter{
    state $counter = 0;
    $counter++;
    if ($counter % 50_000 == 0){
        if ($counter > 1_000_000){
            say STDERR sprintf("%0.3f\n", $counter / 1_000_000);
        }
        else{
            say STDERR $counter;
        }
    }
}
