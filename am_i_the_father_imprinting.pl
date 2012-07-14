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
use lib "$FindBin::Bin/lib";
use Fasta qw/bisulfite_convert/;
use FastqReader;
use Launch;
use DZUtil qw/split_names chext timestamp/;
use Parallel::ForkManager;
my $pm = Parallel::ForkManager->new($opt_parallel);
use Run::BowtieBuild;
use Run::Bowtie;
use MethylCounter;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
unless $opt_output_directory && $opt_reference && $opt_raw && scalar %opt_left_splice;

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
# Handle options

$logger->info("raw file: $opt_raw");
$logger->info("reference A: $opt_reference");
$logger->info("left splice: " . join ',', @opt_left_splice{qw/start end/});
if ($single_sided){
    $logger->info("singled sided (no right splice)");
}
else{
    $logger->info("right splice: " . join ',', @opt_right_splice{qw/start end/});
}
$logger->info("parallel: $opt_parallel");

my $singlecdir_c2t = catfile($opt_output_directory, "single-c-c2t");
my $singlecdir_g2a = catfile($opt_output_directory, "single-c-g2a");
my $freqfile_c2t   = catfile($opt_output_directory, "freq-c2t.txt");
my $freqfile_g2a   = catfile($opt_output_directory, "freq-g2a.txt");
my $windowdir_c2t  = catfile($opt_output_directory, "windows-c2t");
my $windowdir_g2a  = catfile($opt_output_directory, "windows-g2a");

for ( $singlecdir_c2t, $singlecdir_g2a, $windowdir_c2t, $windowdir_g2a){
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
# bs-bowtie-build indices

$logger->info("bs-bowtie-building");

my ($bsrc_reference_c2t, $index_c2t) = bowtie_build(file => $opt_reference, bs => 'c2t', noref => 1);
my ($bsrc_reference_g2a, $index_g2a) = bowtie_build(file => $opt_reference, bs => 'g2a', noref => 1);

#######################################################################
# convert reads

$logger->info("converting fastq->fasta with bs conversion");

my $rawfas = "$opt_raw.fasta";
my $raw_c2t = "$opt_raw.fasta.c2t";
my $raw_g2a = "$opt_raw.fasta.g2a";

FastqReader::fastq_to_fasta(undef, $opt_raw, $rawfas) unless (-f $rawfas && (-s $rawfas) > (-s $opt_raw) / 2);
FastqReader::fastq_to_fasta('c2t', $opt_raw, $raw_c2t) unless (-f $raw_c2t && (-s $raw_c2t) > (-s $opt_raw) / 2);
FastqReader::fastq_to_fasta('g2a', $opt_raw, $raw_g2a) unless (-f $raw_g2a && (-s $raw_g2a) > (-s $opt_raw) / 2);

#######################################################################
# Run Bowtie

$logger->info("running $opt_raw against eco a and b bowtie");
my $basename_c2t = "$basename-c2t";
my $basename_g2a = "$basename-g2a";

my $left_basename_c2t = "$basename_c2t.left";
my $left_basename_g2a = "$basename_g2a.left";

# my $right_basename_c2t = "$basename_c2t.right";
# my $right_basename_g2a = "$basename_g2a.right";

my $left_bowtie_c2t = "$left_basename_c2t.0.bowtie";
my $left_bowtie_g2a = "$left_basename_g2a.0.bowtie";

# my $right_bowtie_c2t = "$right_basename_c2t.0.bowtie";
# my $right_bowtie_g2a = "$right_basename_g2a.0.bowtie";

for my $pair (
    [$raw_c2t, $left_bowtie_c2t, $index_c2t], 
    [$raw_g2a, $left_bowtie_g2a, $index_g2a],
) {
    $logger->info("bowtie $pair->[0] against $pair->[2] into $pair->[1]");
    
    if ($pm->start == 0){
        if (-f $pair->[1] && -s $pair->[1]){
            $logger->info("seems to be done already, skipping bowtie-ing for $pair->[1]");
        }
        else{
            my (undef, undef , undef, undef, @loglines) = bowtie(
                verbose    => 1,
                '-1'       => $pair->[0], 
                output     => $pair->[1],
                index      => $pair->[2],
                splice     => [$opt_left_splice{start}, $opt_left_splice{end}],
                format     => 'fasta',
                norc       => 1,
                base       => 1,
                maxhits    => $opt_max_hits,
                mismatches => $opt_bowtie_mismatches,
            );
            $logger->info($_) for @loglines;

        }
        $pm->finish;
    }
}
$pm->wait_all_children;
# if (!$single_sided){ }

#######################################################################
# parse bowtie

$logger->info("parse bowtie -> eland");

my $left_eland_c2t = "$left_basename_c2t.1.eland";
my $left_eland_g2a = "$left_basename_g2a.1.eland";

# my $right_eland_a = "$right_basename_c2t.1.eland";
# my $right_eland_b = "$right_basename_g2a.1.eland";

if ($pm->start == 0){
    launch("perl -S parse_bowtie.pl $whole_flag -u $rawfas -s @opt_left_splice{qw/start end/} $left_bowtie_c2t -o ??",
        expected => $left_eland_c2t, force => $opt_force);
    $pm->finish;
}
if ($pm->start == 0){
    launch("perl -S parse_bowtie.pl $whole_flag -u $rawfas -s @opt_left_splice{qw/start end/} $left_bowtie_g2a -o ??",
        expected => $left_eland_g2a, force => $opt_force);
    $pm->finish;
}
# if (!$single_sided){
#     if ($pm->start == 0){
#         launch("perl -S parse_bowtie.pl $whole_flag -u $rawfas -s @opt_right_splice{qw/start end/} $right_bowtie_a -o ??",
#             expected => $right_eland_a, force => $opt_force);
#         $pm->finish;
#     }
#     if ($pm->start == 0){
#         launch("perl -S parse_bowtie.pl $whole_flag -u $rawfas -s @opt_right_splice{qw/start end/} $right_bowtie_b -o ??",
#             expected => $right_eland_b, force => $opt_force);
#         $pm->finish;
#     }
# 
# }
$pm->wait_all_children;

#######################################################################
# Split on mismatches

$logger->info("split_on_mismatch.pl");
my $left_eland_filtered_c2t = "$left_basename_c2t.2.elfiltered";
my $left_eland_filtered_g2a = "$left_basename_g2a.2.elfiltered";
# my $right_eland_filtered_c2t = "$right_basename_c2t.2.elfiltered";
# my $right_eland_filtered_g2a = "$right_basename_g2a.2.elfiltered";

launch("perl -S split_on_mismatches_2.pl -a $left_eland_c2t -b $left_eland_g2a -oa $left_eland_filtered_c2t -ob $left_eland_filtered_g2a",
    expected => [ $left_eland_filtered_c2t, $left_eland_filtered_g2a]);

# if (! $single_sided){
#     launch("perl -S split_on_mismatches_2.pl -a $right_eland_a -b $right_eland_b -oa $right_eland_filtered_c2t -ob $right_eland_filtered_g2a",
#         expected => [ $right_eland_filtered_c2t, $right_eland_filtered_g2a]);
# }

#######################################################################
# Unionize left and right

my $eland_union_c2t = "$basename_c2t.3.union";
my $eland_union_g2a = "$basename_g2a.3.union";

if ($single_sided){
    $eland_union_c2t = $left_eland_filtered_c2t;
    $eland_union_g2a = $left_eland_filtered_g2a;
}
# else {
#     launch("perl -S eland_unionize.pl -l $left_eland_filtered_c2t -r $right_eland_filtered_c2t -o $eland_union_c2t", expected => $eland_union_c2t);
#     launch("perl -S eland_unionize.pl -l $left_eland_filtered_g2a -r $right_eland_filtered_g2a -o $eland_union_g2a", expected => $eland_union_g2a);
# }


#######################################################################
# Count and calc ratios and stuff

my $left_ratio  = "$basename.ratio.left.txt";
my $right_ratio = "$basename.ratio.right.txt";
my $ratio       = "$basename.ratio.txt";

# if double sided, calculate the ratios for the left/right individually as well.
# if (!$single_sided){
#     launch("perl -S split_ratio.pl -r $opt_reference_a -o $left_ratio -ea c2t -eb g2a " 
#         . " -a $left_eland_filtered_c2t -b $left_eland_filtered_g2a -m $opt_bowtie_mismatches",
#         expected => "$left_ratio");
# 
#     launch("perl -S split_ratio.pl -r $opt_reference_a -o $right_ratio -ea c2t -eb g2a " 
#         . " -a $right_eland_filtered_c2t -b $right_eland_filtered_g2a -m $opt_bowtie_mismatches",
#         expected => "$right_ratio");
# }

launch("perl -S split_ratio.pl -r $opt_reference -o $ratio -ea c2t -eb g2a -a $eland_union_c2t -b $eland_union_g2a -m $opt_bowtie_mismatches",
    expected => "$ratio");

if (! $opt_no_fracmeth){
    #######################################################################
    # correlate

    my $gff_c2t = "$basename_c2t.4.gff";
    my $gff_g2a = "$basename_g2a.4.gff";

    # make sure reads map together

    if ($pm->start == 0){
        launch("perl -S correlateSingleEnds.pl -rnd 1 -e $eland_union_c2t -f $opt_reference -o $gff_c2t", expected => $gff_c2t);
        $pm->finish;
    }
    if ($pm->start == 0){
        launch("perl -S correlateSingleEnds.pl -rnd 1 -e $eland_union_g2a -f $opt_reference -o $gff_g2a", expected => $gff_g2a);
        $pm->finish;
    }
    $pm->wait_all_children;

    #######################################################################
    # count methyl

    MethylCounter::batch(
        dinucleotide      => 0,
        genome            => $opt_reference,
        correlation       => $gff_c2t,
        parallel          => $pm,
        verbose           => 1,
        bstype            => 'c2t',
        prefix            => catfile($singlecdir_c2t, "$basename_base-c2t-%s"),
        parallel          => $pm,
    );
    MethylCounter::batch(
        dinucleotide      => 0,
        genome            => $opt_reference,
        correlation       => $gff_g2a,
        parallel          => $pm,
        verbose           => 1,
        bstype            => 'g2a',
        prefix            => catfile($singlecdir_g2a, "$basename_base-g2a-%s"),
        parallel          => $pm,
    );
    $pm->wait_all_children;

    for my $singlec (catfile($singlecdir_c2t, "*gff")) {
        $pm->start and next;
        my $window = catfile($windowdir_c2t, basename($singlec, '.gff') . 'w50.gff');

        launch("perl -S window_by_fixed.pl -w $opt_window_size --reference $opt_reference --output ?? --no-skip $singlec", 
            expected => $window);

        $pm->finish; 
    }
    $pm->wait_all_children;

}

=head1 NAME

bastard_imprinting.pl - BaStard = BiSulfite parent imprinting.

=head1 SYNOPSIS

 bastard_imprinting.pl -r raw.fastq -ea Col -b Ler -ra genome-a.fasta -rb genome-b.fasta -mh 10 -m 2 -s 1 50 -o outdir -b basename

=head1 OPTIONS

=over

=item  -r <file> | --raw <file>

FastQ format reads

=for Euclid
    file.type:        readable

=item  -f <fasta> | --reference <fasta>

Genome reference for A.

=for Euclid 
    fasta.type: readable

=item  -ls <start> <end> | --left-splice <start> <end>

=for Euclid
    start.type:     int
    end.type:     int

=item  -rs <start> <end> | --right-splice <start> <end>

=for Euclid
    start.type:     int
    end.type:     int

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

=item  --whole

=item --no-fracmeth

=item  -w <width> | --window-size <width>

=for Euclid
    width.default:     50

=back

=cut

