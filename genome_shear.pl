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

my $num_reads   = 10_000;
my $read_length = 100;
my $out_file;

my $result = GetOptions (
    "num-reads|n=i"   => \$num_reads,
    "read-length|l=i" => \$read_length,
    "out-file|o=s"    => \$out_file,
);

unless (defined $out_file && @ARGV == 1 && -f $ARGV[0]){
    say "usage: $0 [-n #numreads] [-l read_length] [-o out.fastq] input";
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
    my $rc   = rand() > .50;
    say "\@$seq:$start:$stop:" . ($rc ? '-' : '+');
    say $fr->get($seq, $start, $stop, base => 1, rc => $rc);
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

