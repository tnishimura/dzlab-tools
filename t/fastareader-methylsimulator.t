#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More qw(no_plan);
use Test::Exception;

use FastaReader::MethylSimulator;

my $fasta = "t/data/test.fasta";

# no methylation
{
    my $fr = FastaReader::MethylSimulator->new(
        file     => $fasta,
        slurp    => 1,
        methrate => 0,
        norc     => 0,
        readlen  => 100,
    );

    for (1 .. 100) {
        my ($seq, $start, $stop, $rc, $read, $bsread) = @{$fr->get_read()};
        is($fr->get($seq, $start, $stop, rc => $rc), $read, "nonmeth read");
    }

    ok(0 == keys %{$fr->bsrecord()}, "no bs recorded");
}

# methylation
{
    my $fr = FastaReader::MethylSimulator->new(
        file     => $fasta,
        slurp    => 1,
        methrate => .1,
        norc     => 0,
        readlen  => 100,
    );

    my @bsreads;
    for (1 .. 10000) {
        #my ($seq, $start, $stop, $rc, $read, $bsread) = @{$fr->get_read()};
        push @bsreads, $fr->get_read();
    }

    my $expected_meth = find_methylation(\@bsreads);
    say Dumper $expected_meth;

    is_deeply($fr->bsrecord(), $expected_meth, "methylation");
}

sub find_methylation{
    my $bsreads = shift;

    my %results;

    for my $bsread (@$bsreads) {
        my ($seq, $start, $end, $rc, $read_str, $bsread_str) = @$bsread;
        my @read = split //, $read_str;
        my @bsread = split //, $bsread_str;

        for my $abspos ($start .. $end){
            my $relpos = $abspos - $start;
            if ($rc){
                $abspos = $end - $relpos;
            }
            
            if ($read[$relpos] eq 'C' and $bsread[$relpos] eq 'C'){
                $results{$seq}{$abspos} += $rc == 0 ? 1 : -1;
            }
        }
    }

    return \%results;
}


