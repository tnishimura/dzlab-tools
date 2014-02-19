#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;

sub mean_sd{
    my $file = shift;
    my $p = GFF::Parser->new(file => $file, normalize => 0);

    # http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
    my $n = 0;
    my $mean = 0;
    my $m2 = 0;

    while (defined(my $gff = $p->next)){
        my $score = $gff->score;
        next if ! defined $score;

        $n += 1;
        my $delta = $score - $mean;
        $mean += $delta / $n;
        $m2   += $delta * ($score - $mean);
    }

    my $variance = $m2 / ($n - 1);

    return ($mean, sqrt($variance));
}

sub zscore_normalize{
    my $file = shift;
    my ($mean, $sd) = mean_sd($file);
    say STDERR "mean: $mean, sd: $sd";

    my $p = GFF::Parser->new(file => $file, normalize => 0);
    while (defined(my $gff = $p->next)){
        my $score = $gff->score;
        if (defined $score){
            $gff->score(($score - $mean)/$sd);
        }
        say $gff;
    }
}

my $file = shift or pod2usage(-verbose => 2, -noperldoc => 1);

zscore_normalize($file);

=head1 gff-zscore.pl 

Calculate the Z score of the scores column ((row_score - score_mean)/score_stddev)

 gff-zscore.pl input.gff > output.gff

=cut

