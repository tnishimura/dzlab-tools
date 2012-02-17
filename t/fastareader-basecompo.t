#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More qw(no_plan);
use Test::Exception;

BEGIN{
    use_ok "FastaReader::BaseComposition";
}

my $file = "t/data/test.fasta";
{
    my @chr1parts = split //, "CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCAT";
    my $fr = FastaReader->new(file => $file, slurp => 0);
    my $iter = FastaReader::BaseComposition::make_iterator($fr, 'chr1', 3);
    my $motif_len = 3;

    my @results;
    for my $i (0 .. scalar @chr1parts - $motif_len) {
        my $x = join "", @chr1parts[$i .. $i + $motif_len - 1];
        my $y = join "", @{$iter->()};
        push @results, $x eq $y;
    }
    ok(! (grep { ! $_ } @results), "make_iterator()" );
}

{
    my $fr = FastaReader->new(file => $file, slurp => 0);
    my $mc1 = base_composition($file, 1);
    my $mc2 = base_composition($file, 2);
    my $mc3 = base_composition($file, 3);
    my $meth = methylation_context($file);
    my $bc = $fr->base_composition();

    for my $seq ($fr->sequence_list()) {
        for my $x (qw/A C G T/) {
            is_deeply($mc1->{$seq}{$x}, $bc->{$seq}{$x}, "single base composition $seq $x");
        }

        for my $x (keys %{$mc2->{$seq}}){
            is_deeply($mc2->{$seq}{$x}, scalar(()=$fr->get($seq)=~m/(?=$x)/g), "dual base compo $seq $x");
        }

        for my $x (keys %{$mc3->{$seq}}){
            is_deeply($mc3->{$seq}{$x}, scalar(()=$fr->get($seq)=~m/(?=$x)/g), "triple base compo $seq $x");
        }
    }
}

FastaReader::BaseComposition::report($file, 3);
