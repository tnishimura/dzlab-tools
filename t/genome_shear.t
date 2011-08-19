#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use FastaReader;
use Carp;

use Test::More qw(no_plan);

my $fasta = "t/test.fasta";
my $fr = FastaReader->new(file => $fasta, slurp => 1);
my $numreads = 1000;
my $readlen = 100;

open my $shearfh, "./genome_shear.pl -n $numreads -l $readlen -o - $fasta |";

my $counter = 1;
while (! eof $shearfh){
    my @lines = map { scalar <$shearfh> } (1 .. 4);
    chomp @lines;

    # @chrm:35:134:+ 71,78,108
    if ($lines[0] =~ /@(\w+):(\d+):(\d+):(\+|\-)(?:\s([\w,]+))?/){
        my ($seqname, $start, $end, $strand, $meth) = ($1, $2, $3, $4, $5);

        my $read_c2t = c2t($lines[1]);
        my $genome_c2t = c2t($fr->get($seqname, $start, $end, rc => $strand eq '-', base => 1));
        is($genome_c2t, $read_c2t, "read check $counter");

        my $numc = $lines[1] =~ tr/C//;
        if ($meth){
            my @methpos = split /,/, $meth;
            my $nummeth = scalar @methpos;
            is($numc, $nummeth, "read C count check $counter");
            for my $pos (@methpos) {
                is('C', $fr->get($seqname, $pos, $pos, rc => $strand eq '-', base => 1), "pos check $counter");
            }
        }
        else{
            is($numc, 0, "unmeth read $counter has no C's");

        }
    }
    else {
        croak("malformated output from genome_shear.pl?\n" . join "\n", @lines );
    }
    $counter++;
}

close $shearfh;

sub c2t{
    (my $c2t = $_[0]) =~ s/C/T/gi;
    return $c2t;
}

ok(1);
