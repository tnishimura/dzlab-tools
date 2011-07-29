#!/usr/bin/env perl
use 5.10.0;
use strict;
use warnings FATAL => "all";
use Data::Dumper;
use feature 'say';
use autodie;
use Pod::Usage;
use Getopt::Long;
use List::Util qw/sum/;
use FindBin;
use lib "$FindBin::Bin/lib";
use FastaReader;

END {close STDOUT}

# read in genome, bisulte treat at random rate,
# shear into pieces and spit out fastq

my $result = GetOptions (
    "num-reads|n=i"   => \(my $num_reads = 10_000),
    "read-length|l=i" => \(my $read_length = 100),
    "out-file|o=s"    => \(my $out_file),
    "tag|t=s"         => \(my $tag),
    "no-rc|f"         => \(my $norc),
    "reverse-coordinate" => \(my $reverse_coordinates)
);

if ($reverse_coordinates && ! $norc){
    say "--reverse-coordinate only really makes sense with --no-rc, quitting";
    exit 1;
}

unless (defined $out_file && @ARGV == 1 && -f $ARGV[0]){
    say "usage: $0 [--reverse-coordinates] [--no-rc] [-n #numreads] [-l read_length] [-o out.fastq] input";
    say "  --reverse-coordinates flips the coordinates and strand symbol in output. use for rc genomes";
    say "  --no-rc means don't shear the rc strand. use for rc genomes";
    exit 1;
}

if ($out_file ne '-'){
    open my $tmp, '>', $out_file;
    select $tmp;
}

my $fr = FastaReader->new(file => $ARGV[0], slurp => 1);
my %seqlen = $fr->sequence_lengths();
my $multargs = [(),()];

while (my ($seq,$len) = each %seqlen) {
    push @{$multargs->[0]}, $seq;
    push @{$multargs->[1]}, $len;
}

my $quality = "#" x $read_length;

for (1 .. $num_reads){
    my $seq   = rmultinomial(@$multargs);
    my $start = int(rand($seqlen{$seq}-$read_length))+1;
    my $stop  = $start + $read_length -1;
    my $get_rc   = $norc ? 0 : rand() > .50;
    my $strand_label = ($norc && $reverse_coordinates || !$reverse_coordinates && $get_rc) ? '-' : '+';
    
    my ($real_start, $real_stop) = ($start, $stop);
    if ($reverse_coordinates){
        $real_start = ($seqlen{$seq} - $stop  + 1); # flip start/stop
        $real_stop  = ($seqlen{$seq} - $start + 1);
    }

    say "\@$seq:$real_start:$real_stop:$strand_label";
    say $fr->get($seq, $start, $stop, base => 1, rc => $get_rc);
    say "+";
    say $quality;
} 

sub rmultinomial{
    die "rmultinomial requires at least one probability-item pairs" if @_ < 2;
    my ($items, $probabilities) = @_;
    if (ref $items ne 'ARRAY' || ref $probabilities ne 'ARRAY' || @$items != @$probabilities){
        die "usage: rmultinomial [qw/A C G T/], [.25, .25, .25, .25]";
    }

    my $count = @$probabilities;
    my $norm= sum @$probabilities;
    @$probabilities = map { $_ / $norm } @$probabilities;

    if ($count == 1){
        return $items->[0];
    }

    my $unif = rand(1);
    
    my $cummulative_probability = 0;

    for my $i (0 .. $count - 1){

        $cummulative_probability += $probabilities->[$i];
        if ($unif < $cummulative_probability){
            return $items->[$i];
        }
    }
    
    # there's a chance we end up here b/c of floating point issues
    return $items->[-1];
}

