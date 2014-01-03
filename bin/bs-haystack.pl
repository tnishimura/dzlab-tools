#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use File::Basename qw/basename/;
use Getopt::Long qw/:config no_ignore_case/;
use List::MoreUtils qw/any/;
use Pod::Usage;
use IO::All;

use FindBin;
use lib "$FindBin::Bin/../lib";


use FastqReader;
use DZUtil qw/fastq_read_length/;
use Digest::MD5::Util;
use Run::BowtieBuild;
use Run::Bowtie2;

my $result = GetOptions (
    "output-directory|d=s" => \(my $output_directory),
    "reference-file|r=s"   => \(my $reference_file),
    "tdna-file|t=s"        => \(my $tdna_file),
    "reads-file|R=s"       => \(my $reads_file),
);

usage("malformed arguments") if (!$result);
usage("reference, tdna, flank, reads file all need to exist") 
if any { ! defined $_ || ! -f $_ } ($reference_file, $tdna_file, $reads_file);
usage("need --output-directory (-d)") if ! $output_directory;

#######################################################################
# setup output directory and files

my $outdir = io($output_directory)->mkpath();

my $reads_vs_tdna              = $outdir->catfile(basename($reads_file, qw{.fq .fastq}) . "_vs_" .  basename($tdna_file, '.fasta') . ".sam");
my $reads_vs_tdna_mapping      = $reads_vs_tdna . ".mapping";
my $reads_vs_reference         = $outdir->catfile(basename($reads_file, qw{.fq .fastq}) . "_vs_" .  basename($reference_file, '.fasta') . ".sam");
my $reads_vs_reference_mapping = $reads_vs_reference . ".mapping";

my $log_file                = $outdir->catfile("log.txt");
my $reads_bsrc_file         = $outdir->catfile(basename($reads_file) . ".c2t");
my $reads_bsrc_file_mapping = $reads_bsrc_file . ".mapping-to-tdna";

#######################################################################

LOG("additional bowtie2 options: @ARGV");

# convert reads c2t
launch($reads_bsrc_file, "perl -S fastq2rcfasta.pl --c2t $reads_file >  $reads_bsrc_file ");

# bs-bowtie reference 
my ($reference_bsrc_file) = bowtie_build( file => $reference_file, bs => 'c2t', rc => 1, version => 2);

# and scaffold
my ($tdna_bsrc_file) = bowtie_build( file => $tdna_file, bs => 'c2t', rc => 1, version => 2);

# align reads to tdna 
launch($reads_vs_tdna, "bowtie2 --local -f --norc -x $tdna_bsrc_file -U $reads_bsrc_file -S $reads_vs_tdna @ARGV 2>> $log_file");
launch($reads_vs_tdna_mapping,"sam-only-mapped.pl < $reads_vs_tdna > $reads_vs_tdna_mapping");

# grab aligning read's ID's
my %read_id_mapping_to_tdna = map { 
    my @cols = split /\t/, $_;
    @cols >= 10 ? ($cols[0] => 1) : ()
} io($reads_vs_tdna_mapping)->chomp->getlines();

# grab reads from 
LOG("grabbing read ID's from $reads_bsrc_file > $reads_bsrc_file_mapping");
{
    my $in = io($reads_bsrc_file);
    my $out = io($reads_bsrc_file_mapping);

    my $counter = 0;
    while (defined(my $header = $in->getline())){
        LOG($counter) if ++$counter % 2500000 == 0;
        my $body = $in->getline();

        if ($header =~ />([^\s]+)/ and exists $read_id_mapping_to_tdna{$1} ){
            $out->print($header);
            $out->print($body);
        }
    }
};

launch($reads_vs_reference, "bowtie2 --local -f --norc -x $reference_bsrc_file -U $reads_bsrc_file_mapping -S $reads_vs_reference @ARGV 2>> $log_file");
launch($reads_vs_reference_mapping, "sam-only-mapped.pl < $reads_vs_reference > $reads_vs_reference_mapping");

#######################################################################
# aux function

sub launch{
    my $expected = shift;
    my $cmd = shift;
    if (!io($expected)->exists){
        LOG("Running $cmd");
    warn $cmd;
        system($cmd);
    }
    else{
        LOG("$expected already exists, not running $cmd");
    }
}

sub LOG { 
    my $msg = shift;
    io("=")->println("LOG: $msg\n");
    io($log_file)->appendln("LOG: $msg\n");
}
sub usage { 
    pod2usage(scalar(@_) ? (-msg => shift()) : (), -verbose => 2, -noperldoc => 1); 
}


