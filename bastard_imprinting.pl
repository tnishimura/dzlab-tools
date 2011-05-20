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
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/lib";
use Fasta qw/bisulfite_convert/;
use Launch;
use DZUtil qw/split_names chext timestamp/;
use Parallel::ForkManager;
my $pm = Parallel::ForkManager->new($opt_parallel);


pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
unless $opt_output_directory && $opt_reference_a && $opt_reference_b && $opt_raw && $opt_ecotype_a && $opt_ecotype_b && scalar %opt_left_splice;

my $logname = $opt_output_directory . "-" . timestamp() . ".log.txt";

my $conf=qq/
    log4perl.logger          = DEBUG, Print
    log4perl.logger.PipeLine = DEBUG, File

    log4perl.appender.Print                          = Log::Log4perl::Appender::Screen
    log4perl.appender.Print.layout                   = PatternLayout
    log4perl.appender.Print.layout.ConversionPattern = %d{HH:mm:ss} %p> %m%n

    log4perl.appender.File                          = Log::Log4perl::Appender::File
    log4perl.appender.File.filename                 = $logname
    log4perl.appender.File.layout                   = PatternLayout
    log4perl.appender.File.layout.ConversionPattern = %d{HH:mm:ss} %p> (%L) %M - %m%n
/;
Log::Log4perl::init( \$conf );

my $logger = get_logger($opt_debug ? "" : "PipeLine");

my $single_sided = ! scalar %opt_right_splice;

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

my $singlecdir_a = catfile($opt_output_directory, "single-c-$opt_ecotype_a");
my $singlecdir_b = catfile($opt_output_directory, "single-c-$opt_ecotype_b");
my $windowdir_a = catfile($opt_output_directory, "windows-$opt_ecotype_a");
my $windowdir_b = catfile($opt_output_directory, "windows-$opt_ecotype_b");

mkdir $opt_output_directory;
mkdir $singlecdir_a;
mkdir $singlecdir_b;
mkdir $windowdir_a;
mkdir $windowdir_b;
if (! -d $opt_output_directory){ $logger->logdie("can't create $opt_output_directory"); }
if (! -d $singlecdir_a){ $logger->logdie("can't create $singlecdir_a"); }
if (! -d $singlecdir_b){ $logger->logdie("can't create $singlecdir_b"); }
if (! -d $windowdir_a){ $logger->logdie("can't create $windowdir_a"); }
if (! -d $windowdir_b){ $logger->logdie("can't create $windowdir_b"); }

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

my $bsrc_reference_a = $opt_reference_a . ".bsrc";
my $bsrc_reference_b = $opt_reference_b . ".bsrc";

launch("perl -S fasta_bsrc.pl $opt_reference_a > $bsrc_reference_a", expected => $bsrc_reference_a, force => $opt_force >= 2 );
launch("perl -S fasta_bsrc.pl $opt_reference_b > $bsrc_reference_b", expected => $bsrc_reference_b, force => $opt_force >= 2 );

for my $ref ($bsrc_reference_a,$bsrc_reference_b) {
    launch("bowtie-build $ref $ref", expected => ("$ref.1.ebwt"), force => $opt_force >= 2);
}

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

if ($pm->start == 0){
    launch("bowtie $bsrc_reference_a -f -B 1 -v $opt_bowtie_mismatches --norc --best -5 $left_trim5 -3 $left_trim3 $rawc2t $left_bowtie_a",
        expected => $left_bowtie_a, force => $opt_force);
    $pm->finish;
}
if ($pm->start == 0){
    launch("bowtie $bsrc_reference_b -f -B 1 -v $opt_bowtie_mismatches --norc --best -5 $left_trim5 -3 $left_trim3 $rawc2t $left_bowtie_b",
        expected => $left_bowtie_b, force => $opt_force);
    $pm->finish;
}

if (!$single_sided){
    if ($pm->start == 0){
        launch("bowtie $bsrc_reference_a -f -B 1 -v $opt_bowtie_mismatches --norc --best -5 $right_trim5 -3 $right_trim3 $rawc2t $right_bowtie_a",
            expected => $right_bowtie_a, force => $opt_force);
        $pm->finish;
    }
    if ($pm->start == 0){
        launch("bowtie $bsrc_reference_b -f -B 1 -v $opt_bowtie_mismatches --norc --best -5 $right_trim5 -3 $right_trim3 $rawc2t $right_bowtie_b",
            expected => $right_bowtie_b, force => $opt_force);
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
    launch("perl -S parse_bowtie.pl -w -u $rawfas -s @opt_left_splice{qw/start end/} $left_bowtie_a -o $left_eland_a",
        expected => $left_eland_a, force => $opt_force);
    $pm->finish;
}
if ($pm->start == 0){
    launch("perl -S parse_bowtie.pl -w -u $rawfas -s @opt_left_splice{qw/start end/} $left_bowtie_b -o $left_eland_b",
        expected => $left_eland_b, force => $opt_force);
    $pm->finish;
}
if (!$single_sided){
    if ($pm->start == 0){
        launch("perl -S parse_bowtie.pl -w -u $rawfas -s @opt_right_splice{qw/start end/} $right_bowtie_a -o $right_eland_a",
            expected => $right_eland_a, force => $opt_force);
        $pm->finish;
    }
    if ($pm->start == 0){
        launch("perl -S parse_bowtie.pl -w -u $rawfas -s @opt_right_splice{qw/start end/} $right_bowtie_b -o $right_eland_b",
            expected => $right_eland_b, force => $opt_force);
        $pm->finish;
    }

}
$pm->wait_all_children;

#######################################################################
# Split on mismatches

$logger->info("split_on_mismatch.pl");
my $left_eland_filtered_a = "$left_basename_a.2.elfiltered";
my $left_eland_filtered_b = "$left_basename_b.2.elfiltered";
my $right_eland_filtered_a = "$right_basename_a.2.elfiltered";
my $right_eland_filtered_b = "$right_basename_b.2.elfiltered";

launch("perl -S split_on_mismatches_2.pl -a $left_eland_a -b $left_eland_b -oa $left_eland_filtered_a -ob $left_eland_filtered_b",
    expected => [ $left_eland_filtered_a, $left_eland_filtered_b]);

if (! $single_sided){
    launch("perl -S split_on_mismatches_2.pl -a $right_eland_a -b $right_eland_b -oa $right_eland_filtered_a -ob $right_eland_filtered_b",
        expected => [ $right_eland_filtered_a, $right_eland_filtered_b]);
}

#######################################################################
# Unionize left and right

my $eland_union_a = "$basename_a.3.union";
my $eland_union_b = "$basename_b.3.union";

if ($single_sided){
    $eland_union_a = $left_eland_filtered_a;
    $eland_union_b = $left_eland_filtered_b;
}
else {
    launch("eland_unionize.pl -l $left_eland_filtered_a -r $right_eland_filtered_a -o $eland_union_a", expected => $eland_union_a);
    launch("eland_unionize.pl -l $left_eland_filtered_b -r $right_eland_filtered_b -o $eland_union_b", expected => $eland_union_b);
}


#######################################################################
# Count and calc ratios and stuff

# if double sided, calculate the ratios for the left/right individually as well.
if (!$single_sided){
    launch("perl -S split_ratio.pl -o $basename.ratio.left.txt -ea $opt_ecotype_a -eb $opt_ecotype_b " 
        . " -a $left_eland_filtered_a -b $left_eland_filtered_b -m $opt_bowtie_mismatches",
        expected => "$basename.ratio.left.txt");

    launch("perl -S split_ratio.pl -o $basename.ratio.right.txt -ea $opt_ecotype_a -eb $opt_ecotype_b " 
        . " -a $right_eland_filtered_a -b $right_eland_filtered_b -m $opt_bowtie_mismatches",
        expected => "$basename.ratio.right.txt");
}

launch("perl -S split_ratio.pl -o $basename.ratio.txt -ea $opt_ecotype_a -eb $opt_ecotype_b -a $eland_union_a -b $eland_union_b -m $opt_bowtie_mismatches",
    expected => "$basename.ratio.txt");


if (! $opt_no_fracmeth){
    #######################################################################
    # correlate

    my $gff_a = "$basename_a.4.gff";
    my $gff_b = "$basename_b.4.gff";

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

    launch("perl -S countMethylation_batch.pl -g $gff_a -s $singlecdir_a -w $windowdir_a -r $opt_reference_a -b $basename_base-vs-$opt_ecotype_a -t $opt_parallel");
    launch("perl -S countMethylation_batch.pl -g $gff_b -s $singlecdir_b -w $windowdir_b -r $opt_reference_b -b $basename_base-vs-$opt_ecotype_b -t $opt_parallel");

    launch("collect-freqs.pl -o $basename.freq $opt_output_directory");
}

#######################################################################
# split_on_mismatches WITH coord checking. 

if ($opt_check_coord_ratio){
    # Split on mismatches
    $logger->info("split_on_mismatch.pl");
    my $cc_left_eland_filtered_a = "$left_basename_a.2.check_coord.elfiltered";
    my $cc_left_eland_filtered_b = "$left_basename_b.2.check_coord.elfiltered";
    my $cc_right_eland_filtered_a = "$right_basename_a.2.check_coord.elfiltered";
    my $cc_right_eland_filtered_b = "$right_basename_b.2.check_coord.elfiltered";

    launch("perl -S split_on_mismatches_2.pl --check-coord -a $left_eland_a -b $left_eland_b -oa $cc_left_eland_filtered_a -ob $cc_left_eland_filtered_b",
        expected => [ $cc_left_eland_filtered_a, $cc_left_eland_filtered_b]);

    if (! $single_sided){
        launch("perl -S split_on_mismatches_2.pl --check-coord -a $right_eland_a -b $right_eland_b -oa $cc_right_eland_filtered_a -ob $cc_right_eland_filtered_b",
            expected => [ $cc_right_eland_filtered_a, $cc_right_eland_filtered_b]);
    }

    # Unionize left and right
    my $cc_eland_union_a = "$basename_a.3.check_coord.union";
    my $cc_eland_union_b = "$basename_b.3.check_coord.union";

    if ($single_sided){
        $cc_eland_union_a = $cc_left_eland_filtered_a;
        $cc_eland_union_b = $cc_left_eland_filtered_b;
    }
    else {
        launch("eland_unionize.pl -l $cc_left_eland_filtered_a -r $cc_right_eland_filtered_a -o $cc_eland_union_a", expected => $cc_eland_union_a);
        launch("eland_unionize.pl -l $cc_left_eland_filtered_b -r $cc_right_eland_filtered_b -o $cc_eland_union_b", expected => $cc_eland_union_b);
    }

    # Count and calc ratios and stuff

    # if double sided, calculate the ratios for the left/right individually as well.
    if (!$single_sided){
        launch("perl -S split_ratio.pl -o $basename.ratio.check_coord.left.txt -ea $opt_ecotype_a -eb $opt_ecotype_b " 
            . " -a $cc_left_eland_filtered_a -b $cc_left_eland_filtered_b -m $opt_bowtie_mismatches",
            expected => "$basename.ratio.check_coord.left.txt");

        launch("perl -S split_ratio.pl -o $basename.ratio.check_coord.right.txt -ea $opt_ecotype_a -eb $opt_ecotype_b " 
            . " -a $cc_right_eland_filtered_a -b $cc_right_eland_filtered_b -m $opt_bowtie_mismatches",
            expected => "$basename.ratio.check_coord.right.txt");
    }

    launch("perl -S split_ratio.pl -o $basename.ratio.check_coord.txt -ea $opt_ecotype_a -eb $opt_ecotype_b -a $cc_eland_union_a -b $cc_eland_union_b -m $opt_bowtie_mismatches",
        expected => "$basename.ratio.check_coord.txt");

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

=item --debug

=item --help | -h

=item  --force <level>

Level of forcefulness in doing jobs.  1 = Redo all run-specifics.  2 = Redo bowtie-build as well.

=for Euclid
    level.default:     0
    level.type:        int, level >= 0

=item  --merge 

=item  -cc | --check-coord-ratio

Enabling this option will create a second set of filtered eland and resulting
ratio files with coordinate checking on split_on_mismatches_2.pl enabled.  You
probably don't need this. 

=item --no-fracmeth

=back

=cut

