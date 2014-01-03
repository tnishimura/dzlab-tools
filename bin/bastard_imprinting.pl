#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use List::Util qw/sum/;
use Log::Log4perl qw/get_logger/;
use File::Spec::Functions;
use File::Basename;
use File::Path;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Fasta qw/bisulfite_convert/;
use Launch;
use DZUtil qw/split_names chext timestamp/;
use Parallel::ForkManager;
my $pm = Parallel::ForkManager->new($opt_parallel);

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
unless $opt_output_directory && $opt_reference_a && $opt_reference_b && $opt_raw && $opt_ecotype_a && $opt_ecotype_b && scalar %opt_left_splice;

if (! -d $opt_output_directory) { mkpath ( $opt_output_directory ,  {verbose => 1} ); }

my $logname = catfile($opt_output_directory, $opt_basename) . "-" . timestamp() . ".log.txt";
my $bowtie_logname = catfile($opt_output_directory, $opt_basename). "-" . timestamp() . ".log-bowtie.txt";

my $conf=qq/
    log4perl.logger          = DEBUG, Print
    log4perl.logger.PipeLine = DEBUG, File

    log4perl.appender.Print                          = Log::Log4perl::Appender::Screen
    log4perl.appender.Print.layout                   = PatternLayout
    log4perl.appender.Print.layout.ConversionPattern = %d{HH:mm:ss} %.1p> %m%n

    log4perl.appender.File                          = Log::Log4perl::Appender::File
    log4perl.appender.File.filename                 = $logname
    log4perl.appender.File.layout                   = PatternLayout
    log4perl.appender.File.layout.ConversionPattern = %d{HH:mm:ss} %.1p> - %m%n
/;
Log::Log4perl::init( \$conf );

my $logger = get_logger($opt_debug ? "" : "PipeLine");

my $single_sided = ! scalar %opt_right_splice;
my $whole_flag = $opt_whole ? '-w' : '';

#################################################################################

my ($left_trim5, $left_trim3)   = ($opt_left_splice{start} - 1,  $opt_read_length - $opt_left_splice{end});
my ($right_trim5, $right_trim3); 
if (!$single_sided){
    ($right_trim5, $right_trim3) = ($opt_right_splice{start} - 1, $opt_read_length - $opt_right_splice{end});
}

#######################################################################
# Handle options


$logger->info("raw file: $opt_raw");
$logger->info("reference A: $opt_reference_a");
$logger->info("reference B: $opt_reference_b");
$logger->info("ecotype A: $opt_ecotype_a");
$logger->info("ecotype B: $opt_ecotype_b");
$logger->info("left splice: " . join ',', @opt_left_splice{qw/start end/});
if ($single_sided){
    $logger->info("singled sided (no right splice)");
}
else{
    $logger->info("right splice: " . join ',', @opt_right_splice{qw/start end/});
}
$logger->info("parallel: $opt_parallel");

my $nocc = $opt_no_cc ? ".no_coord_check" : "";
my $nocc_flag = $opt_no_cc ? "" : "--check-coord";

my $singlecdir_a = catfile($opt_output_directory, "single-c-$opt_ecotype_a$nocc");
my $singlecdir_b = catfile($opt_output_directory, "single-c-$opt_ecotype_b$nocc");
my $windowdir_a  = catfile($opt_output_directory, "windows-$opt_ecotype_a$nocc");
my $windowdir_b  = catfile($opt_output_directory, "windows-$opt_ecotype_b$nocc");

for ( $singlecdir_a, $singlecdir_b, $windowdir_a, $windowdir_b){
    mkdir $_;
    if (! -d $_) {$logger->logdie("can't create $_");}
}

my $basename = $opt_basename;
if (! defined $basename){
    $basename = $opt_raw;
    $basename =~ s/\.\w+$//;
    $basename = basename $basename;
}
my $basename_base = $basename;
$basename = catfile($opt_output_directory,$basename);

$logger->info("basename: $basename");

#######################################################################
# BSRC genomes, build bowtie indices

my $bsrc_reference_a = $opt_reference_a . ".c2t";
my $bsrc_reference_b = $opt_reference_b . ".c2t";

#launch("perl -S fasta_bsrc.pl $opt_reference_a > $bsrc_reference_a", expected => $bsrc_reference_a, force => $opt_force >= 2 );
#launch("perl -S fasta_bsrc.pl $opt_reference_b > $bsrc_reference_b", expected => $bsrc_reference_b, force => $opt_force >= 2 );

launch("perl -S bs-bowtie-build -c2t $opt_reference_a");
launch("perl -S bs-bowtie-build -c2t $opt_reference_b");

#for my $ref ($bsrc_reference_a,$bsrc_reference_b) {
#    launch("bowtie-build $ref $ref", expected => ("$ref.1.ebwt"), force => $opt_force >= 2);
#}

# raw: fastq -> fasta,  c2t
$logger->info("converting fastq->fasta, then converting");
my $rawfas = "$opt_raw.fasta";
my $rawc2t = "$opt_raw.fasta.c2t";
launch("perl -S fq_all2std.pl fq2fa $opt_raw > $rawfas", expected => $rawfas, force => $opt_force>=2);
launch("perl -S convert.pl c2t $rawfas > $rawc2t", expected => $rawc2t, force => $opt_force>=2);

#######################################################################
# Run Bowtie

# run bowtie of rawc2t vs output a and b
$logger->info("running $opt_raw against eco a and b bowtie");
my $basename_a = "$basename-vs-$opt_ecotype_a";
my $basename_b = "$basename-vs-$opt_ecotype_b";

my $left_basename_a = "$basename_a.left";
my $left_basename_b = "$basename_b.left";

my $right_basename_a = "$basename_a.right";
my $right_basename_b = "$basename_b.right";

my $left_bowtie_a = "$left_basename_a.0.bowtie";
my $left_bowtie_b = "$left_basename_b.0.bowtie";

my $right_bowtie_a = "$right_basename_a.0.bowtie";
my $right_bowtie_b = "$right_basename_b.0.bowtie";

my $max_hits_arg = $opt_max_hits ? " --strata -k 1 -m $opt_max_hits " : "";

if ($pm->start == 0){
    launch("bowtie $bsrc_reference_a $max_hits_arg -f -B 1 -v $opt_bowtie_mismatches --norc --best -5 $left_trim5 -3 $left_trim3 $rawc2t ??",
        also => $bowtie_logname, expected => $left_bowtie_a, force => $opt_force);
    $pm->finish;
}
if ($pm->start == 0){
    launch("bowtie $bsrc_reference_b $max_hits_arg -f -B 1 -v $opt_bowtie_mismatches --norc --best -5 $left_trim5 -3 $left_trim3 $rawc2t ??",
        also => $bowtie_logname, expected => $left_bowtie_b, force => $opt_force);
    $pm->finish;
}

if (!$single_sided){
    if ($pm->start == 0){
        launch("bowtie $bsrc_reference_a $max_hits_arg -f -B 1 -v $opt_bowtie_mismatches --norc --best -5 $right_trim5 -3 $right_trim3 $rawc2t ??",
            also => $bowtie_logname, expected => $right_bowtie_a, force => $opt_force, );
        $pm->finish;
    }
    if ($pm->start == 0){
        launch("bowtie $bsrc_reference_b $max_hits_arg -f -B 1 -v $opt_bowtie_mismatches --norc --best -5 $right_trim5 -3 $right_trim3 $rawc2t ??",
            also => $bowtie_logname, expected => $right_bowtie_b, force => $opt_force, );
        $pm->finish;
    }
}

$pm->wait_all_children;

#######################################################################
# parse bowtie

$logger->info("parse bowtie -> eland");

my $left_eland_a = "$left_basename_a.1.eland";
my $left_eland_b = "$left_basename_b.1.eland";

my $right_eland_a = "$right_basename_a.1.eland";
my $right_eland_b = "$right_basename_b.1.eland";


if ($pm->start == 0){
    launch("perl -S parse_bowtie.pl $whole_flag -u $rawfas -s @opt_left_splice{qw/start end/} $left_bowtie_a -o ??",
        expected => $left_eland_a, force => $opt_force);
    $pm->finish;
}
if ($pm->start == 0){
    launch("perl -S parse_bowtie.pl $whole_flag -u $rawfas -s @opt_left_splice{qw/start end/} $left_bowtie_b -o ??",
        expected => $left_eland_b, force => $opt_force);
    $pm->finish;
}
if (!$single_sided){
    if ($pm->start == 0){
        launch("perl -S parse_bowtie.pl $whole_flag -u $rawfas -s @opt_right_splice{qw/start end/} $right_bowtie_a -o ??",
            expected => $right_eland_a, force => $opt_force);
        $pm->finish;
    }
    if ($pm->start == 0){
        launch("perl -S parse_bowtie.pl $whole_flag -u $rawfas -s @opt_right_splice{qw/start end/} $right_bowtie_b -o ??",
            expected => $right_eland_b, force => $opt_force);
        $pm->finish;
    }

}
$pm->wait_all_children;

#######################################################################
# Split on mismatches


$logger->info("split_on_mismatch.pl");
my $left_eland_filtered_a = "$left_basename_a.2.elfiltered$nocc";
my $left_eland_filtered_b = "$left_basename_b.2.elfiltered$nocc";
my $right_eland_filtered_a = "$right_basename_a.2.elfiltered$nocc";
my $right_eland_filtered_b = "$right_basename_b.2.elfiltered$nocc";

launch("perl -S split_on_mismatches_2.pl $nocc_flag -a $left_eland_a -b $left_eland_b -oa $left_eland_filtered_a -ob $left_eland_filtered_b",
    expected => [ $left_eland_filtered_a, $left_eland_filtered_b]);

if (! $single_sided){
    launch("perl -S split_on_mismatches_2.pl $nocc_flag -a $right_eland_a -b $right_eland_b -oa $right_eland_filtered_a -ob $right_eland_filtered_b",
        expected => [ $right_eland_filtered_a, $right_eland_filtered_b]);
}

#######################################################################
# Unionize left and right

my $eland_union_a = "$basename_a.3.union$nocc";
my $eland_union_b = "$basename_b.3.union$nocc";

if ($single_sided){
    $eland_union_a = $left_eland_filtered_a;
    $eland_union_b = $left_eland_filtered_b;
}
else {
    launch("perl -S eland_unionize.pl -l $left_eland_filtered_a -r $right_eland_filtered_a -o $eland_union_a", expected => $eland_union_a);
    launch("perl -S eland_unionize.pl -l $left_eland_filtered_b -r $right_eland_filtered_b -o $eland_union_b", expected => $eland_union_b);
}


#######################################################################
# Count and calc ratios and stuff

my $left_ratio  = "$basename.ratio.left.txt$nocc";
my $right_ratio = "$basename.ratio.right.txt$nocc";
my $ratio       = "$basename.ratio.txt$nocc";

# if double sided, calculate the ratios for the left/right individually as well.
if (!$single_sided){
    launch("perl -S split_ratio.pl -r $opt_reference_a -o $left_ratio -ea $opt_ecotype_a -eb $opt_ecotype_b " 
        . " -a $left_eland_filtered_a -b $left_eland_filtered_b -m $opt_bowtie_mismatches",
        expected => "$left_ratio");

    launch("perl -S split_ratio.pl -r $opt_reference_a -o $right_ratio -ea $opt_ecotype_a -eb $opt_ecotype_b " 
        . " -a $right_eland_filtered_a -b $right_eland_filtered_b -m $opt_bowtie_mismatches",
        expected => "$right_ratio");
}

launch("perl -S split_ratio.pl -r $opt_reference_a -o $ratio -ea $opt_ecotype_a -eb $opt_ecotype_b -a $eland_union_a -b $eland_union_b -m $opt_bowtie_mismatches",
    expected => "$ratio");

if (! $opt_no_fracmeth){
    #######################################################################
    # correlate

    my $gff_a = "$basename_a.4.gff$nocc";
    my $gff_b = "$basename_b.4.gff$nocc";

    # make sure reads map together

    if ($pm->start == 0){
        launch("perl -S correlateSingleEnds.pl -rnd 1 -e $eland_union_a -f $opt_reference_a -o $gff_a", expected => $gff_a);
        $pm->finish;
    }
    if ($pm->start == 0){
        launch("perl -S correlateSingleEnds.pl -rnd 1 -e $eland_union_b -f $opt_reference_b -o $gff_b", expected => $gff_b);
        $pm->finish;
    }
    $pm->wait_all_children;

    #######################################################################
    # count methyl

    if ($pm->start == 0){
        launch("perl -S countMethylation_batch.pl -g $gff_a -s $singlecdir_a -w $windowdir_a -r $opt_reference_a -b $basename_base-vs-$opt_ecotype_a$nocc");
        $pm->finish;
    }
    if ($pm->start == 0){
        launch("perl -S countMethylation_batch.pl -g $gff_b -s $singlecdir_b -w $windowdir_b -r $opt_reference_b -b $basename_base-vs-$opt_ecotype_b$nocc");
        $pm->finish;
    }
    $pm->wait_all_children;

    launch("perl -S collect-freqs.pl -o $basename.freq $opt_output_directory");
}

=head1 NAME

bastard_imprinting.pl - BaStard = BiSulfite parent imprinting.

=head1 SYNOPSIS

Usage examples:

 bastard_imprinting.pl -r raw.fastq -ea Col -b Ler -ra genome-a.fasta -rb genome-b.fasta -mh 10 -l 100 -m 2 -s 1 50 -o outdir -b basename

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

=item  -ea <eco> | --ecotype-a <eco>

Ecotype A label.

=item  -eb <eco> | --ecotype-b <eco>

Ecotype B label.

=item  -ls <start> <end> | --left-splice <start> <end>

=for Euclid
    start.type:     int
    end.type:     int

=item  -rs <start> <end> | --right-splice <start> <end>

=for Euclid
    start.type:     int
    end.type:     int

=item  -l <len> | --read-length <len>

Number of bp in each read.

=item  -m <num> | --bowtie-mismatches <num>

Number of mismatches to allow in bowtie.

=for Euclid
    num.default:     2
    num.type:        int, num >= 0 && num <= 3

=item -o <dir> | --output-directory <dir>

Output Directory.

=item  -b <name> | --basename <name>

Prefix for the file names.

=item  --parallel <numthreads>

Number of concurrent jobs to run.  CAREFUL.  Default 0, for no parallelization. 

=for Euclid
    numthreads.default:     0

=item  -mh <hits> | --max-hits <hits>

=item --debug

=item --help | -h

=item  --force <level>

Level of forcefulness in doing jobs.  1 = Redo all run-specifics.  2 = Redo bowtie-build as well.

=for Euclid
    level.default:     0
    level.type:        int, level >= 0

=item  --merge 

=item  --no-cc

Enabling this option will create a second set of filtered eland and resulting
ratio files with coordinate checking on split_on_mismatches_2.pl disabled.  You
probably don't need this. 

=item  -w | --whole

=item --no-fracmeth

=back

=cut

