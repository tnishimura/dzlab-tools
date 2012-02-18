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

TODO: {

    my $fr = FastaReader->new(file => $file, slurp => 0);
    my ($single, $compo) = base_composition($file, 3);
    my $bc = $fr->base_composition();

    #say Dumper $mc3;
    for my $seq ($fr->sequence_list()) {
        for my $x (qw/A C G T/) {
            is_deeply($single->{$seq}{$x}, $bc->{$seq}{$x}, "single base composition $seq $x");
        }

        my @contexts = keys %{$compo->{$seq}};
        for my $x (grep { !($_ ~~ [qw/A C G T CHG CHH/]) } @contexts){
            is_deeply($compo->{$seq}{$x}, scalar(()=$fr->get($seq)=~m/(?=$x)/g), "base compo $seq $x");
        }
    }
}

{
    my $fr = FastaReader->new(file => $file, slurp => 0);
    my ($single, $compo) = base_composition($file, 3, 1);
    #say Dumper $mc3;
    for my $seq ($fr->sequence_list()) {
        my %contexts = (CG => qr/CG/, CHG => qr/C[ACT]G/, CHH => qr/C[ACT][ACT]/);
        while (my ($name,$rx) = each %contexts) {
            my $forward = scalar(()=$fr->get($seq)=~m/(?=$rx)/g);
            my $reverse = scalar(()=$fr->get($seq, undef, undef, rc => 1)=~m/(?=$rx)/g);

            #say "$forward, $reverse";
            is_deeply($compo->{$seq}{$name}, $forward + $reverse, "methyl base composition $seq $name");
        }
    }
}


