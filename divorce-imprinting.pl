#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use List::Util qw/sum/;
use File::Spec::Functions;
use File::Basename;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/lib";
use Launch;
use DZUtil qw/chext timestamp/;
use Parallel::ForkManager;
my $pm = Parallel::ForkManager->new(2);

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
unless $opt_output_directory && $opt_annotation && $opt_reference_a && $opt_reference_b && $opt_raw && $opt_ecotype_a && $opt_ecotype_b && scalar %opt_splice;

my $logname = $opt_output_directory . "-" . timestamp() . ".log.txt";

use Log::Log4perl qw/get_logger/;
my $conf=qq/
    log4perl.logger          = DEBUG, Print
    log4perl.logger.PipeLine = DEBUG, File

    log4perl.appender.Print        = Log::Log4perl::Appender::Screen
    log4perl.appender.Print.layout = PatternLayout
    log4perl.appender.Print.layout.ConversionPattern = %d{HH:mm:ss} %p> (%L) %M - %m%n

    log4perl.appender.File          = Log::Log4perl::Appender::File
    log4perl.appender.File.filename = $logname
    log4perl.appender.File.layout   = PatternLayout
    log4perl.appender.File.syswrite = 1
    log4perl.appender.File.layout.ConversionPattern = %d{HH:mm:ss} %p> (%L) %M - %m%n
/;
Log::Log4perl::init( \$conf );
my $logger = get_logger($opt_debug ? "" : "PipeLine");
#my $logger = get_logger();

#################################################################################

my ($trim5, $trim3) = ($opt_splice{start} - 1, $opt_read_length - $opt_splice{end});

#######################################################################
# Handle options

$logger->info("raw file: $opt_raw");
$logger->info("reference A: $opt_reference_a");
$logger->info("reference B: $opt_reference_b");
$logger->info("ecotype A: $opt_ecotype_a");
$logger->info("ecotype B: $opt_ecotype_b");
$logger->info("splice: " . join ',', @opt_splice{qw/start end/});

mkdir $opt_output_directory;
if (! -d $opt_output_directory){
    $logger->logdie("can't create $opt_output_directory");
}

my $basename = $opt_basename;
if (! defined $basename){
    $basename = $opt_raw;
    $basename =~ s/\.\w+$//;
    $basename = basename $basename;
}
$basename = catfile($opt_output_directory,$basename);

$logger->info("basename: $basename");

#######################################################################
# BSRC genomes, build bowtie indices

for my $ref ($opt_reference_a,$opt_reference_b) {
    launch("bowtie-build $ref $ref", expected => ("$ref.1.ebwt"), force => $opt_force >= 2);
}

# raw: fastq -> fasta,  c2t
$logger->info("converting fastq->fasta, then converting");
my $rawfas = chext($opt_raw,"fasta");

launch("perl -S fq_all2std.pl fq2fa $opt_raw > $rawfas", expected => $rawfas, force => $opt_force>=2);

#######################################################################
# Run Bowtie

# run bowtie of rawc2t vs output a and b
$logger->info("running $opt_raw against eco a and b bowtie");
my $basename_a = "$basename-vs-$opt_ecotype_a";
my $basename_b = "$basename-vs-$opt_ecotype_b";

my $bowtie_a = "$basename_a.0.bowtie";
my $bowtie_b = "$basename_b.0.bowtie";

if ($pm->start == 0){
    launch("bowtie $opt_reference_a -f -B 1 -v $opt_bowtie_mismatches --best -5 $trim5 -3 $trim3 $rawfas $bowtie_a",
        expected => $bowtie_a, force => $opt_force);
    $pm->finish();
}
if ($pm->start == 0){
    launch("bowtie $opt_reference_b -f -B 1 -v $opt_bowtie_mismatches --best -5 $trim5 -3 $trim3 $rawfas $bowtie_b",
        expected => $bowtie_b, force => $opt_force);
    $pm->finish();
}
$pm->wait_all_children;


#######################################################################
# parse bowtie

$logger->info("parse bowtie -> eland");
my $eland_a = "$basename_a.1.eland";
my $eland_b = "$basename_b.1.eland";

if ($pm->start == 0){
    launch("perl -S parse_bowtie.pl -u $rawfas -s @opt_splice{qw/start end/} $bowtie_a -o $eland_a",
        expected => $eland_a, force => $opt_force);
    $pm->finish();
}
if ($pm->start == 0){
    launch("perl -S parse_bowtie.pl -u $rawfas -s @opt_splice{qw/start end/} $bowtie_b -o $eland_b",
        expected => $eland_b, force => $opt_force);
    $pm->finish();
}
$pm->wait_all_children;

#######################################################################
# Sort

$logger->info("sort the eland files");
my $eland_sorted_a = "$basename_a.2.elsorted";
my $eland_sorted_b = "$basename_b.2.elsorted";

if ($pm->start == 0){
    launch("sort -k 1,1 -k 4,4 -S 100M $eland_a -o $eland_sorted_a",expected => $eland_sorted_a, force => $opt_force);
    $pm->finish();
}
if ($pm->start == 0){
    launch("sort -k 1,1 -k 4,4 -S 100M $eland_b -o $eland_sorted_b",expected => $eland_sorted_b, force => $opt_force);
    $pm->finish();
}
$pm->wait_all_children;

#######################################################################
# Split on mismatches

$logger->info("split_on_mismatch.pl");
my $eland_filtered_a = "$basename_a.3.elfiltered";
my $eland_filtered_b = "$basename_b.3.elfiltered";

launch("perl -S split_on_mismatches_2.pl -a $eland_a -b $eland_b -oa $eland_filtered_a -ob $eland_filtered_b",
    expected => [ $eland_filtered_a, $eland_filtered_b]);

#######################################################################
# Parse_eland.pl

$logger->info("Parse eland to gff");

my $gff_a = "$basename_a.4.gff";
my $gff_b = "$basename_b.4.gff";


if ($pm->start == 0){
    launch("perl -S parse_eland.pl -3 $eland_filtered_a -o $gff_a", expected => $gff_a);
    $pm->finish();
}
if ($pm->start == 0){
    launch("perl -S parse_eland.pl -3 $eland_filtered_b -o $gff_b", expected => $gff_b);
    $pm->finish();
}
$pm->wait_all_children;

#######################################################################
# Parse_eland.pl

$logger->info("sort the gff's");

my $gff_sorted_a = "$basename_a.5.sorted.gff";
my $gff_sorted_b = "$basename_b.5.sorted.gff";

if ($pm->start == 0){
    launch("sort -k 1,1 -k 4,4n -k 5,5n -k 7,7 -S 100M $gff_a -o $gff_sorted_a", expected => $gff_sorted_a);
    $pm->finish();
}
if ($pm->start == 0){
    launch("sort -k 1,1 -k 4,4n -k 5,5n -k 7,7 -S 100M $gff_b -o $gff_sorted_b", expected => $gff_sorted_b);
    $pm->finish();
}
$pm->wait_all_children;

#######################################################################
# filter_gff

$logger->info("filter_repeats");

my $gff_filtered_a = "$basename_a.6.filtered.gff";
my $gff_filtered_b = "$basename_b.6.filtered.gff";

my $gff_repeats_a = "$basename_a.6.repeats.gff";
my $gff_repeats_b = "$basename_b.6.repeats.gff";

if ($pm->start == 0){
    launch("perl -S filter_repeats.pl $gff_sorted_a -o $gff_filtered_a 2> $gff_repeats_a", expected => $gff_filtered_a);
    $pm->finish();
}
if ($pm->start == 0){
    launch("perl -S filter_repeats.pl $gff_sorted_b -o $gff_filtered_b 2> $gff_repeats_b", expected => $gff_filtered_b);
    $pm->finish();
}
$pm->wait_all_children;

#######################################################################
# windowing

$logger->info("filter_repeats");

my $w50_a = "$basename_a.7.w50.gff";
my $w50_b = "$basename_b.7.w50.gff";

my $w50_filtered_a = "$basename_a.7.w50-filtered.gff";
my $w50_filtered_b = "$basename_b.7.w50-filtered.gff";

if ($pm->start == 0){
    launch("perl -S window_gff.pl $gff_sorted_a -g $opt_annotation -k -c sum -o $w50_a -r", expected => $w50_a);
    launch("perl -S window_gff.pl $gff_filtered_a -g $opt_annotation -k -c sum -o $w50_filtered_a -r", expected => $w50_filtered_a);
    $pm->finish();
}
if ($pm->start == 0){
    launch("perl -S window_gff.pl $gff_sorted_b -g $opt_annotation -k -c sum -o $w50_b -r", expected => $w50_b);
    launch("perl -S window_gff.pl $gff_filtered_b -g $opt_annotation -k -c sum -o $w50_filtered_b -r", expected => $w50_filtered_b);
    $pm->finish();
}
$pm->wait_all_children;

my $table = "$basename.table.txt";

launch("perl -S divorce_gene_table.pl -a $opt_annotation -f $w50_a $w50_filtered_a $w50_b $w50_filtered_b -o $table", expected => $table);

=head1 NAME

ratio.pl - Your program here

=head1 SYNOPSIS

Usage examples:

 ratio.pl -r raw.fastq -ea Col -b Ler -ra genome-a.fasta -rb genome-b.fasta -l 100 -m 2 -s 1 50 -o outdir -b basename

=head1 REQUIRED ARGUMENTS

=over

=back

=head1 OPTIONS

=over

=item  -r <file> | --raw <file>

FastQ format reads

=for Euclid
    file.type:        readable

=item  -ra <fasta> | --reference-a <fasta>

Genome reference for A.

=for Euclid 
    fasta.type: readable

=item  -rb <fasta> | --reference-b <fasta>

Genome reference for B.

=for Euclid 
    fasta.type: readable

=item  -a <a> | --annotation <a>

GFF annotation file

=for Euclid
    a.type:        readable

=item  -ea <eco> | --ecotype-a <eco>

Ecotype A label.

=item  -eb <eco> | --ecotype-b <eco>

Ecotype B label.

=item  -s <start> <end> | --splice <start> <end>

=for Euclid
    start.type:     int
    end.type:     int

=item  -l <len> | --read-length <len>

Number of bp in each read.

=for Euclid
    len.default: 100

=item  -m <num> | --bowtie-mismatches <num>

Number of mismatches to allow in bowtie.

=for Euclid
    num.default:     2
    num.type:        int, num >= 0 && num <= 3

=item -o <dir> | --output-directory <dir>

Output Directory.

=item  -b <name> | --basename <name>

Prefix for the file names.

=item --help | -h

=item  --force <level>

Level of forcefulness in doing jobs.  1 = Redo all run-specifics.  2 = Redo bowtie-build as well.

=for Euclid
    level.default:     0
    level.type:        int, level >= 0

=item  --debug

=back

=cut

