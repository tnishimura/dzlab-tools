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
    "bs-rate|r=f"     => \(my $bsrate = .1),
    "out-file|o=s"    => \(my $out_file),
    "tag|t=s"         => \(my $tag),
    "no-rc|f"         => \(my $norc),
    "junk|j=i"        => \(my $junk=0),
);

unless (defined $out_file && @ARGV == 1 && -f $ARGV[0]){
    say "usage: $0 [--no-rc] [--bs-rate .1] [-j #junk] [-n #numreads] [-l read_length] [-o out.fastq] input";
    say "  --no-rc means don't shear the rc strand. use for rc genomes";
    exit 1;
}

my $logfile;
if ($out_file ne '-'){
    open my $tmp, '>', $out_file;
    select $tmp;
    $logfile = $out_file . ".methsites.gff";
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
    my $strand_label = $get_rc ? '-' : '+';
    
    my $read = $fr->get($seq, $start, $stop, base => 1, rc => $get_rc);
    my ($bsread, $meth) = bisulfite($read, $start, $stop, $get_rc, $bsrate);

    record_meth($seq, $strand_label, @$meth);

    say "\@$seq:$start:$stop:$strand_label:" , join(":",@$meth); 
    say $bsread;
    say "+";
    say $read;
} 

for (1 .. $junk){
    say "\@junk$_";
    say 'N' x $read_length;
    say "+";
    say 'N' x $read_length;

}

if ($logfile){
    dump_meth($logfile);
}

{
    my %bsrecord; 
    sub record_meth{
        my ($seq, $strand, @meth) = @_;
        for my $pos (@meth) {
            # if forward strand, positive. if rc strand, negative. I'm so crafty.
            if (($strand eq '+' && exists $bsrecord{$seq}{$pos} && $bsrecord{$seq}{$pos} < 0) || 
                ($strand eq '-' && exists $bsrecord{$seq}{$pos} && $bsrecord{$seq}{$pos} > 0)){
                die "methylation on BOTH sides? $seq, $pos";
            }
            $bsrecord{$seq}{$pos} += $strand eq '+' ? 1 : -1;
        }
    }
    sub dump_meth{
        open my $fh, '>', $logfile;
        for my $seq (sort keys %bsrecord) {
            for my $pos (sort {$a <=> $b} keys %{$bsrecord{$seq}}){
                my $score = $bsrecord{$seq}{$pos};
                say $fh join "\t", $seq, qw/. C/, $pos, $pos, abs($score), $score > 0 ? '+' : '-', qw/. ./;
            }
        }
    }
}

sub bisulfite{
    my ($subseq, $start, $end, $getrc, $rate) = @_;
    my @split = split //, $subseq;
    my @meth;
    for my $abspos ($start .. $end){
        my $relpos = $abspos - $start;
        if ($getrc){
            $abspos = $end - $relpos;
        }

        if ($split[$relpos] eq 'C'){
            if( rand() > $rate){
                $split[$relpos] = 'T';
            }
            else{
                push @meth, $abspos;
            }
        }
    }
    @meth = sort { $a <=> $b } @meth;
    return join("", @split), \@meth; 
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

