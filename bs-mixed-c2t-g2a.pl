#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use 5.010_000;
use Carp;
use List::Util qw/sum/;
use Log::Dispatch;
use File::Spec::Functions;
use File::Basename;
use File::Path;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

use FindBin;
use lib "$FindBin::Bin/lib";
use Fasta qw/bisulfite_convert/;
use FastqReader;
use FastqReader::Convert;
use Launch;
use DZUtil qw/split_names chext timestamp/;
use Parallel::ForkManager;
use Run::BowtieBuild;
use Run::Bowtie;
use MethylCounter;
use GFF::Statistics qw/methyl_stats/;

my $pm = Parallel::ForkManager->new($opt_parallel);

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
unless $opt_read_length && $opt_basename && $opt_output_directory && $opt_reference && $opt_raw && scalar %opt_left_splice;

if (! -d $opt_output_directory) { mkpath ( $opt_output_directory ,  {verbose => 1} ); }

my $logname = catfile($opt_output_directory, $opt_basename) . "-" . timestamp() . ".log.txt";
my $bowtie_logname = catfile($opt_output_directory, $opt_basename). "-" . timestamp() . ".log-bowtie.txt";

{
    my $logger = Log::Dispatch->new(
        outputs => [
            [ 'File',   newline => 1, min_level => 'debug', filename => $logname ],
            [ 'Screen', newline => 1, min_level => 'debug' ],
        ],
    );
    sub LOG { $logger->log(level => 'info', message => shift) }
    sub LOGDIE { $logger->log_and_die(level => 'info', message => shift) }
}

my $single_sided = ! scalar %opt_right_splice;
my $whole_flag = $opt_whole ? '-w' : '';

#################################################################################
# Handle options

LOG("raw file: $opt_raw");
LOG("reference A: $opt_reference");
LOG("left splice: " . join ',', @opt_left_splice{qw/start end/});
if ($single_sided){
    LOG("singled sided (no right splice)");
}
else{
    LOG("right splice: " . join ',', @opt_right_splice{qw/start end/});
}
LOG("parallel: $opt_parallel");

my $singlecdir_c2t = catfile($opt_output_directory, "single-c-c2t");
my $singlecdir_g2a = catfile($opt_output_directory, "single-c-g2a");
my $freqfile_c2t   = catfile($opt_output_directory, "freq-c2t.txt");
my $freqfile_g2a   = catfile($opt_output_directory, "freq-g2a.txt");
my $windowdir_c2t  = catfile($opt_output_directory, "windows-c2t");
my $windowdir_g2a  = catfile($opt_output_directory, "windows-g2a");

for ( $singlecdir_c2t, $singlecdir_g2a, $windowdir_c2t, $windowdir_g2a){
    mkdir $_;
    if (! -d $_) {LOGDIE("can't create $_");}
}

my $basename = $opt_basename;
if (! defined $basename){
    $basename = $opt_raw;
    $basename =~ s/\.\w+$//;
    $basename = basename $basename;
}
my $basename_base = $basename;
$basename = catfile($opt_output_directory,$basename);

LOG("basename: $basename");

#######################################################################
# bs-bowtie-build indices

LOG("bs-bowtie-building");

my ($bsrc_reference_c2t, $index_c2t) = bowtie_build(file => $opt_reference, bs => 'c2t', noref => 1);
my ($bsrc_reference_g2a, $index_g2a) = bowtie_build(file => $opt_reference, bs => 'g2a', noref => 1);

#######################################################################
# convert reads

LOG("converting fastq->fasta with bs conversion");

my $rawfas = "$opt_raw.fasta";
my $raw_c2t = "$opt_raw.fasta.c2t";
my $raw_g2a = "$opt_raw.fasta.g2a";

# tends to be very slow, so do it in parallel.
if ($pm->start == 0){
    fastq_convert(in => $opt_raw, out => $rawfas) unless (-f $rawfas && (-s $rawfas) > (-s $opt_raw) / 2);
    $pm->finish();
}
if ($pm->start == 0){
    fastq_convert(methyl => 'c2t', in => $opt_raw, out => $raw_c2t) unless (-f $raw_c2t && (-s $raw_c2t) > (-s $opt_raw) / 2);
    $pm->finish();
}
if ($pm->start == 0){
    fastq_convert(methyl => 'g2a', in => $opt_raw, out => $raw_g2a) unless (-f $raw_g2a && (-s $raw_g2a) > (-s $opt_raw) / 2);
    $pm->finish();
}
$pm->wait_all_children;

#######################################################################
# Run Bowtie

LOG("running $opt_raw against eco a and b bowtie");
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
    LOG("bowtie $pair->[0] against $pair->[2] into $pair->[1]");
    
    if ($pm->start == 0){
        if (-f $pair->[1] && -s $pair->[1]){
            LOG("seems to be done already, skipping bowtie-ing for $pair->[1]");
        }
        else{
            my (undef, undef , undef, undef, @loglines) = bowtie(
                verbose    => 1,
                '-1'       => $pair->[0], 
                output     => $pair->[1],
                index      => $pair->[2],
                readlength => $opt_read_length,
                splice     => [$opt_left_splice{start}, $opt_left_splice{end}],
                format     => 'fasta',
                norc       => 1,
                base       => 1,
                maxhits    => $opt_max_hits,
                mismatches => $opt_bowtie_mismatches,
            );
            LOG($_) for @loglines;

        }
        $pm->finish;
    }
}
$pm->wait_all_children;
# if (!$single_sided){ }

#######################################################################
# parse bowtie

LOG("parse bowtie -> eland");

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

LOG("split_on_mismatch.pl");
my $left_eland_filtered_c2t = "$left_basename_c2t.2.elfiltered";
my $left_eland_filtered_g2a = "$left_basename_g2a.2.elfiltered";

launch("perl -S split_on_mismatches_2.pl -a $left_eland_c2t -b $left_eland_g2a -oa $left_eland_filtered_c2t -ob $left_eland_filtered_g2a",
    expected => [ $left_eland_filtered_c2t, $left_eland_filtered_g2a]);

#######################################################################
# Unionize left and right

my $eland_union_c2t = "$basename_c2t.3.union";
my $eland_union_g2a = "$basename_g2a.3.union";

if ($single_sided){
    $eland_union_c2t = $left_eland_filtered_c2t;
    $eland_union_g2a = $left_eland_filtered_g2a;
}

#######################################################################
# Count and calc ratios and stuff

my $left_ratio  = "$basename.ratio.left.txt";
my $right_ratio = "$basename.ratio.right.txt";
my $ratio       = "$basename.ratio.txt";

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

    #######################################################################
    # methylation stats

    my $methyl_stats_file = "$basename.mstats.txt";
    if (! -f $methyl_stats_file || ! -s $methyl_stats_file){
        # gff-methyl-stats
        my @c2t_singlec_concat = grep { /ALL/ } glob catfile($singlecdir_c2t, "*gff");
        my @g2a_singlec_concat = grep { /ALL/ } glob catfile($singlecdir_g2a, "*gff");

        my (undef, $methyl_stat_output) = methyl_stats( map {
            basename($_,'.gff') => $_
            } @c2t_singlec_concat, @g2a_singlec_concat
        );
        open my $fh, '>', $methyl_stats_file;
        say $fh $methyl_stat_output;
        close $fh;
    }

    #######################################################################
    # window methylation

    for my $singlec (glob catfile($singlecdir_c2t, "*gff")) {
        next if $singlec =~ /ALL/;
        $pm->start and next;
        my $window = catfile($windowdir_c2t, basename($singlec, '.gff') . '.w50.gff');
        launch("perl -S window_by_fixed.pl -w $opt_window_size --reference $opt_reference --output ?? --no-skip $singlec", 
            expected => $window);
        $pm->finish; 
    }
    for my $singlec (glob catfile($singlecdir_g2a, "*gff")) {
        next if $singlec =~ /ALL/;
        $pm->start and next;
        my $window = catfile($windowdir_g2a, basename($singlec, '.gff') . '.w50.gff');
        launch("perl -S window_by_fixed.pl -w $opt_window_size --reference $opt_reference --output ?? --no-skip $singlec", 
            expected => $window);
        $pm->finish; 
    }
    $pm->wait_all_children;

    launch("methyl-coverage.pl -r $opt_reference -p $opt_output_directory/$basename_base.coverage -- --$basename_base $opt_output_directory/single-c-*/*ALL*",
        expected => "$opt_output_directory/$basename_base.coverage.svg",
    );
}

$pm->wait_all_children;

=head1 NAME

bs-mixed-c2t-g2a.pl - like bs-sequel.pl, but for when methylation could be c2t or g2a

=head1 SYNOPSIS

 bs-mixed-c2t-g2a.pl -b basename -rl 100 -r reads.fastq -f reference.fasta -ls 1 50 -m 2 -o outdir -mh 10 --parallel 3

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

=item  --parallel <threads>

=for Euclid
    threads.default:     0
    threads.type:        int, threads >= 0 
    threads.type.error:  <threads> must be positive

=item  -rl <len> | --read-length <len>

=back

=cut

