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
use FastaReader::MethylSimulator;

END {close STDOUT}

# read in genome, bisulte treat at random rate,
# shear into pieces and spit out fastq

my $result = GetOptions (
    "num-reads|n=i"        => \(my $num_reads = 10_000),
    "read-length|l=i"      => \(my $read_length = 100),
    "methylation-rate|r=f" => \(my $methrate),
    "out-file|o=s"         => \(my $out_file),
    "no-rc|f"              => \(my $norc),
    "junk|j=i"             => \(my $junk=0),
);

unless (defined $out_file && @ARGV == 1 && -f $ARGV[0]){
    say "usage: $0 [--no-rc] [--methylation-rate .1] [-j #junk] [-n #numreads] [-l read_length] [-o out.fastq] input";
    say "  --no-rc means don't shear the rc strand. use for rc genomes";
    exit 1;
}

my $logfile;
if ($out_file ne '-'){
    open my $tmp, '>', $out_file;
    select $tmp;
    $logfile = $out_file . ".methsites.gff";
}

my $fr = FastaReader::MethylSimulator->new(
    file        => $ARGV[0],
    slurp       => 1,
    read_length => $read_length,
    methrate    => $methrate,
    norc        => $norc,
);

my $quality = "#" x $read_length;

for (1 .. $num_reads){
    my ($seq, $start, $stop, $rc, $read, $bsread, $meth) = @{$fr->get_read()};
    my $strand_label = $rc ? '-' : '+';

    if (defined $methrate){
        say "\@$seq:$start:$stop:$strand_label:" , join(":",@$meth);
        say $bsread;
        say "+";
        say $read;
    }
    else{
        say "\@$seq:$start:$stop:$strand_label";
        say $read;
        say "+";
        say $read;
    }
} 

for (1 .. $junk){
    say "\@junk$_";
    say 'N' x $read_length;
    say "+";
    say 'N' x $read_length;
}

if ($logfile){
    $fr->dump($logfile);
}
