#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Pod::Usage;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/lib";
use BowtieParser;
use BigArray;
use FastaReader;
END {close STDOUT}
$| = 1;

my $result = GetOptions (
    "reference|r=s"       => \my $reference,
    "window-size|w=i" => \(my $window_size = 50),
    "no-skip|k"           => \(my $noskip),
    "base|b=i"            => \(my $base = 1),
);

if (!$result || !$reference){
    say "usage: bowtie_window.pl -r genome --no-skip --base 1 --window-size 50 bowtiefile";
    exit 1;
}

my $bowtie_reader = BowtieParser->new(file => \*ARGV);
my $fasta_reader = FastaReader->new(file => $reference, slurp => 0);

my %counters = map { 
    uc($_) => BigArray->new(base => $base, size => $fasta_reader->get_length($_))
} $fasta_reader->sequence_list();

my $counter = 0;
while (defined(my $bowtie = $bowtie_reader->next())){
    my (undef, undef, $chr, $pos, $read) = @$bowtie;
    $counters{uc $chr}->increment_range($pos, $pos + length($read) - 1);
    say STDERR $counter if (++$counter % 10000 == 0);
}

for my $seq (sort $fasta_reader->sequence_list()) {
    my $ba = $counters{uc $seq};
    my $start = 1;
    my $max = $fasta_reader->get_length($seq);
    #say STDERR "$seq from $start to $max";
    while ($start <= $max){
        if ($window_size == 1){
            my $value = $ba->{pdl}->at($start - $base); 
            if ($value > 0 || $noskip){
                say join "\t", $seq, qw/. ./, $start, $start, $value, qw/. . ./;
            }
            ++$start;
        }
        else{
            my $end = $start + $window_size - 1;
            if ($end > $max){
                $end = $max;
            }
            my $pdl = $ba->get_range($start, $end);
            my $count = $pdl->max();
            if ($count > 0 || $noskip){
                say join "\t", $seq, qw/. ./, $start, $end, $count, qw/. . ./;
            }

            $start += $window_size;
        }
    }
}


=head1 NAME

bowtie_window.pl - return windows, with scores showing how many reads overlap that window

=head1 SYNOPSIS

Usage examples:

 bowtie_window.pl [options]...

=head1 OPTIONS

=over

=item --help | -h

=back

=cut

