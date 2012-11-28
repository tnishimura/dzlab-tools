#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.014_000;
use Data::Dumper;
use autodie;

use FindBin;
use lib "$FindBin::Bin/lib";
use BowtieParser;

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "feature|f=s" => \(my $feature),
    "max-entropy|m=i" => \(my $max_entropy = .17),
);
usage() if !$result || !$feature;

my $p = BowtieParser->new(file => \*ARGV);

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

my %stats;
my %position = (
    A => 1,
    C => 2,
    G => 3,
    T => 4,
);

while (defined(my $alignment = $p->next(1))){
    my ($readid, $strand, $seqid, $start, $read, $quality, $mystery, $mismatches) = @$alignment;
    for my $mm (@$mismatches) {
        my ($abs_coord, $base_in_ref, $base_in_read) = @$mm;

        if ($base_in_read ne 'N'){
            my $record = $stats{$seqid}{$abs_coord};
            if ($record){
                if ($record->[0] eq $base_in_ref){
                    $record->[$position{$base_in_read}]++;
                }
                else{
                    die "Bowtie error";
                }
            }
            else{
                $stats{$seqid}{$abs_coord} = [$base_in_ref, 0, 0, 0, 0];
                $stats{$seqid}{$abs_coord}[$position{$base_in_read}]++;
            }

            #say join "\t", $seqid, $abs_coord, $base_in_ref, $base_in_read;
        }
    }
    counter();
}
for my $seq (keys %stats) {
    for my $pos (keys %{$stats{$seq}}) {
        my ($base_in_ref, $A, $C, $G, $T) = @{$stats{$seq}{$pos}};
        # say join "\t", $seq, $pos, $base_in_ref, entropy($A, $C, $G, $T), $A, $C, $G, $T;
        my $entropy = entropy($A, $C, $G, $T);
        if ($entropy < $max_entropy){
            # score is -log($entropy) so that higher scores = less entropy
            say join "\t", lc($seq), $base_in_ref, $feature, $pos, $pos, -log($entropy), qw/. ./, "A=$A;C=$C;G=$G;T=$T;entropy=$entropy";
        }
    }
}

use List::Util qw/sum/;

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
