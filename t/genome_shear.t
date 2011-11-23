#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use FastaReader;
use Carp;
use Test::More tests => 6;
use List::MoreUtils qw/all/;

my $fasta = "t/test.fasta";
my $fr = FastaReader->new(file => $fasta, slurp => 1);
my $numreads = 1000;
my $readlen = 100;

my $uneven_fastq;
my $uneven_read_length;
my $malformed_header;
my $read_genome_mismatch;
my $wrong_meth;
my $non_c_meth;
my $counter = 1;

open my $shearfh, "./genome_shear.pl -r .1 -n $numreads -l $readlen -o - $fasta |";
OUTER:
while (! eof $shearfh){
    my @lines = map { scalar <$shearfh> } (1 .. 4);
    chomp @lines;
    if (! all {defined} @lines){
        $uneven_fastq = 1;
        last OUTER;
    }
    if (length $lines[1] != length $lines[3]){
        $uneven_read_length = 1;
        last OUTER;
    }

    # @chr3:163702:163801:+:163780
    # TGTGTGTGTTATTTAAAATTAATTAATATTAATAAATTATAGTATATATAAAAATTTTGATATTTGAAGTTTGGTTGACTAGATAATATTTAATAGGTTT
    # +
    # CGTGCGTGTTACTCAAAATTAATTAATACCAATAAATTATAGTATATATAAAAATTTCGACACCTGAAGTTTGGTTGACCAGACAATATTCAATAGGTTT

    if ($lines[0] =~ /@(\w+):(\d+):(\d+):(\+|\-)(?:([\d:]+))?/){
        my ($seqname, $start, $end, $strand, $meth) = ($1, $2, $3, $4, $5);

        # say $lines[0];
        # say join ",", $seqname, $start, $end, $strand, $meth // 'nometh';

        my $read_c2t = c2t($lines[1]);
        my $genome_c2t = c2t($fr->get($seqname, $start, $end, rc => $strand eq '-', base => 1));

        if ($read_c2t ne $genome_c2t){
            say "read doesn't match genome";
            say join "\n", @lines;
            $read_genome_mismatch = 1;
            last OUTER;
        }

        my $numc = $lines[1] =~ tr/C//;
        if ($meth){
            $meth =~ s/^://; # strip off leading :
            my @methpos = split /:/, $meth;
            my $nummeth = scalar @methpos;
            if ($numc != $nummeth){
                say "c count wrong";
                say join "\n", @lines;
                $wrong_meth = 1;
                last OUTER;
            }
            for my $pos (@methpos) {
                if('C' ne $fr->get($seqname, $pos, $pos, rc => $strand eq '-', base => 1)){
                    say "non-c meth position";
                    say join "\n", @lines;
                    $non_c_meth = 1;
                    last OUTER;
                }
            }
        }
        else{
            say "C's exists when there aren't supposed to be any";
            say join "\n", @lines;
            $wrong_meth = 1;
            last OUTER;
        }
    }
    else {
        say "malformed header";
        say join "\n", @lines;
        $malformed_header = 1;
        last OUTER;
    }
    $counter++;
}

close $shearfh;

ok(!$uneven_fastq, "fastq has line count multiple of 4") ;
ok(!$uneven_read_length, "fastq has uneven read length");
ok(!$malformed_header, "headers ok");
ok(!$read_genome_mismatch, "read match the genome");
ok(!$wrong_meth, "meth count ok");
ok(!$non_c_meth, "meth only at c site");

sub c2t{
    (my $c2t = $_[0]) =~ s/C/T/gi;
    return $c2t;
}
