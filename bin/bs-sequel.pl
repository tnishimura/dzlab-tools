#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use File::Spec::Functions;
use File::Path;
use File::Basename;
use FindBin;
use Parallel::ForkManager;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

use File::Copy;
use lib "$FindBin::Bin/../lib";

use DZUtil qw/downsample mfor timestamp split_names fastq_read_length single_c_concat/;
use Launch;
use GFF::Split;
use MethylCounter;

my $pm = Parallel::ForkManager->new($opt_parallel);

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
unless (
    $opt_left_read && $opt_reference && $opt_base_name 
);

# handle some default optios here

if (! defined $opt_out_directory){
    $opt_out_directory = $opt_base_name;
}

if (! -d $opt_out_directory) { mkpath ( $opt_out_directory ,  {verbose => 1} ); }

my $logname = catfile($opt_out_directory, $opt_base_name) . "-" . timestamp() . ".log.txt";
my $bowtie_logname = catfile($opt_out_directory, $opt_base_name). "-" . timestamp() . ".log-bowtie.txt";

use Log::Log4perl qw/get_logger/;
my $conf=qq/
    log4perl.logger          = DEBUG, Print
    log4perl.logger.PipeLine = DEBUG, File

    log4perl.appender.Print        = Log::Log4perl::Appender::Screen
    log4perl.appender.Print.layout = PatternLayout
    log4perl.appender.Print.layout.ConversionPattern = %d{HH:mm:ss} %.1p> %m%n

    log4perl.appender.File          = Log::Log4perl::Appender::File
    log4perl.appender.File.filename = $logname
    log4perl.appender.File.layout   = PatternLayout
    log4perl.appender.File.syswrite = 1
    log4perl.appender.File.layout.ConversionPattern = %d{HH:mm:ss} %.1p> (%L) %m%n
/;
Log::Log4perl::init( \$conf );
my $logger = get_logger("PipeLine");

my $dry = defined $opt_dry;
my $opt_single_ends = ! $opt_right_read;

#######################################################################
# groups (chromosomes)

my @groups        = (); 	# sequences (chr1, chr2, ...)

# Get all chromosomes, pseudo-chromosomes, groups, etc, in the fasta reference file
# Discards all information after first blank character in fasta header
unless (@groups) {
    open my $REFERENCE, '<', $opt_reference or croak "Can't open $opt_reference: $!";

    while (<$REFERENCE>) {
        if (m/>/) {
            tr/A-Z/a-z/;
            m/>([^\s]+)/ && push @groups, $1;
        }
    }
    close $REFERENCE or carp "Can't close $opt_reference: $!";
}

#######################################################################
# Library size

if ($opt_single_ends && !defined $opt_library_size) {
    $opt_library_size = 0;
}
elsif ($opt_single_ends && $opt_library_size != 0){
    $logger->logdie("if single-ends, -k/--library-size needs to be omitted or 0.");
}
elsif (! $opt_single_ends && !defined $opt_library_size){
    $logger->logdie("if paired ends, -k needs to be defined");
}

#######################################################################
# downsample

my $read_size = fastq_read_length($opt_left_read);
if (defined($opt_right_read) && $read_size != fastq_read_length($opt_right_read)){
    $logger->logdie("$opt_right_read and $opt_left_read not same read lengths?");
}

if (defined $opt_downsample){
    $opt_left_read = downsample $opt_left_read, $opt_downsample, 4, dirname $opt_left_read;
    if (defined $opt_right_read){
        $opt_right_read = downsample $opt_right_read, $opt_downsample, 4, dirname $opt_right_read;
    }
}

#######################################################################
# splice argument sanitizing

if (! %opt_left_splice){
    $logger->logwarn("Left splice -ls was not defined, using full length.");
    $opt_left_splice{start} = 1;
    $opt_left_splice{end} = $read_size;
}

my $do_right = %opt_right_splice;

if (! $do_right && !$opt_single_ends){
    $logger->logdie("You can't specify do paired ends without -rs!");
}

# copy to array

my @left_splice  = @opt_left_splice{qw/start end/};
my @right_splice = $do_right ? (@opt_right_splice{qw/start end/}) : ();

# die if out of bounds

if ($left_splice[0] < 1 || $left_splice[1] > $read_size || $left_splice[1] < $left_splice[0]){
    die "left splice out of bounds";
}

if ($do_right && ($right_splice[0] < 1 || $right_splice[1] > $read_size || $right_splice[1] < $right_splice[0])){
    die "left splice out of bounds";
}

#######################################################################
# Log options

$logger->info("These are the options you ran with:");
$logger->info("--reference $opt_reference");
$logger->info("--left-read $opt_left_read");
if (defined $opt_right_read){
    $logger->info("--right-read $opt_right_read");
}
else{
    $logger->info("--right-read <empty>");
}
$logger->info("--left-splice " . join " ", @opt_left_splice{qw/start end/});
if (%opt_right_splice){
    $logger->info("--right-splice " . join " ", @opt_right_splice{qw/start end/});
}
else{
    $logger->info("--right-splice <empty>");
}
$logger->info("--base-name $opt_base_name");
$logger->info("--out-directory $opt_out_directory");
$logger->info("--library-size " . ($opt_library_size // "undef"));
$logger->info("--organism $opt_organism");
$logger->info("--window-size $opt_window_size");
$logger->info("--single-ends $opt_single_ends");
$logger->info("--max-hits $opt_max_hits");
$logger->info("--mismatches $opt_mismatches");
$logger->info("--random-assign $opt_random_assign");
$logger->info("--trust-dash-2 $opt_trust_dash_2");
$logger->info("--di-nuc-freqs $opt_di_nuc_freqs");
$logger->info("--batch $opt_batch");
$logger->info("--no-windowing : " . ($opt_no_windowing ? "<enable>" : "<disable>"));
$logger->info("--parallel $opt_parallel");
$logger->info("--dry : " . ($opt_dry ? "<enable>" : "<disable>"));


#######################################################################
# Directory creation

my $windows_dir  = catfile ($opt_out_directory, 'windows');
my $single_c_dir = catfile ($opt_out_directory, 'single-c');

# check whether output directory should be created, exists, should overwritten
# root dir created at beginning
if (! -d $windows_dir) { mkpath ( $windows_dir,  {verbose => 1} ); }
if (! -d $single_c_dir){ mkpath ( $single_c_dir, {verbose => 1} ); }

my $basename_left  = catfile($opt_out_directory, basename($opt_left_read));

my $basename_right = $do_right && $opt_right_read ? 
                     catfile($opt_out_directory, basename($opt_right_read)) :
                     $do_right && ! $opt_right_read ? 
                     catfile($opt_out_directory, basename($opt_left_read)) :
                     "";

my $basename       = catfile($opt_out_directory, $opt_base_name);

#######################################################################
# Contexts

my @contexts;		# CG, CHG, CHH, etc.
unless (@contexts) {
    if ($opt_di_nuc_freqs) {@contexts = qw(CA CC CG CT)}
    else {@contexts = qw(CG CHG CHH)}
}

my $eland_left_post  = "${basename_left}_$left_splice[0]-$left_splice[1].eland3.post";
my $eland_right_post = $do_right && ! $opt_single_ends ? 
                       "${basename_right}_$right_splice[0]-$right_splice[1].eland3.post" : 
                       $do_right && $opt_single_ends ? 
                       "${basename_left}_$right_splice[0]-$right_splice[1].eland3.post" : 
                       "";

my $int_left  = $opt_bowtie_file ? ("--intermediate $eland_left_post") : '';
my $int_right = $opt_bowtie_file ? ("--intermediate $eland_right_post") : '';

$int_left  =~ s/eland3\.post$/bowtie/;
$int_right =~ s/eland3\.post$/bowtie/;

#######################################################################
# convert genomes & build bowtie indices
# do this here as well as in bs-bowtie, so there's no race condition

launch("perl -S bs-bowtie-build -c2t $opt_reference");
unless ($opt_single_ends) {
    launch("perl -S bs-bowtie-build -g2a $opt_reference");
}


#######################################################################
# bowtie

# align with bowtie
if ($pm->start == 0){
    launch("perl -S bs-bowtie -r $opt_left_read -f $opt_reference -s @left_splice -n $opt_mismatches -mh $opt_max_hits $int_left -o ??",
        expected => $eland_left_post, dryrun => $dry, also => $bowtie_logname);
    $pm->finish;
}

if ($pm->start == 0){
    if ($do_right){
        if ($opt_single_ends) {
            launch("perl -S bs-bowtie -r $opt_left_read -f $opt_reference -s @right_splice -n $opt_mismatches -mh $opt_max_hits $int_right -o ??",
                expected => $eland_right_post, dryrun => $dry, also => $bowtie_logname);
        }
        else {
            launch("perl -S bs-bowtie -g2a -r $opt_right_read -f $opt_reference -s @right_splice -n $opt_mismatches -mh $opt_max_hits $int_right -o ??",
                expected => $eland_right_post, dryrun => $dry, also => $bowtie_logname);
        }
    }
    $pm->finish;
}

$pm->wait_all_children;

#######################################################################
# Correlate

my $base_gff = "$basename.gff";
my $base_log = "$basename.align_stats.log";

if ($do_right){
    launch("perl -S correlatePairedEnds.pl -l $eland_left_post -r $eland_right_post -ref $opt_reference -o ?? -t 0 -d $opt_library_size -s $read_size -2 $opt_trust_dash_2 -1 $opt_single_ends -m $opt_max_hits -a $opt_random_assign", expected =>  $base_gff, dryrun => $dry);
    launch("perl -S collect_align_stats.pl -1 $eland_left_post -2 $eland_right_post -c $base_gff -t $opt_organism -b $opt_batch > ??", expected =>  $base_log, dryrun => $dry);
}
else{
    launch("perl -S correlateSingleEnds.pl -e $eland_left_post --reference $opt_reference -o ?? -m $opt_max_hits --random-assign $opt_random_assign", expected =>  $base_gff, dryrun => $dry);
    launch("perl -S collect_align_stats.pl -1 $eland_left_post -c $base_gff -t $opt_organism -b $opt_batch > ??", expected =>  $base_log, dryrun => $dry);
}

#######################################################################
# Split correlate

$logger->info("Splitting $base_gff by group");

my %base_gff_split;    

if (!$dry){
    %base_gff_split = GFF::Split::split_sequence($base_gff,@groups);
    $logger->info("result of context split of $base_gff: \n" . join "\n", values %base_gff_split);
}
else {
    $logger->info("DRY: GFF::Split::split_sequence($base_gff,@groups);");
}


#######################################################################
# discountMethylation.pl

if ($opt_new_cm){
    $logger->info("starting --new-cm");
    while (my ($seq,$base_gff_part) = each %base_gff_split) {
        $pm->start and next;
        my $single_c_prefix = catfile($single_c_dir, $opt_base_name) . ".$seq";
        my %single = map { $_ => "$single_c_prefix.single-c.$_.gff" } @contexts;
        my $freqfile = "$single_c_prefix.single-c.freq";

        $logger->info("started --new-cm on $seq (creating $single_c_prefix.* and $freqfile)");

        if (4 != grep { -f && -s } $freqfile, values %single){
            my $mc = MethylCounter->new(
                dinucleotide => $opt_di_nuc_freqs,
                genome       => $opt_reference,
                correlation  => $base_gff_part,
                verbose      => $opt_verbose,
            );

            $mc->process();
            $mc->output_single_c( %single );
            $mc->print_freq($freqfile);
        }

        unless ($opt_no_windowing){
            for my $sc (values %single) {
                my $windows_base = basename($sc);
                if ($windows_base !~ s/\.single-c/.w$opt_window_size/){
                    $logger->logdie("$sc- naming screwed up?");
                }
                my $window = catfile($windows_dir, $windows_base);

                launch("perl -S window_by_fixed.pl -w $opt_window_size --reference $opt_reference --output ?? --no-skip $sc", 
                    expected => $window, dryrun => $dry);
            }
        }

        $logger->info("DONE!");
        $pm->finish; 
    }
    $pm->wait_all_children;
    MethylCounter::combine($opt_di_nuc_freqs, "$basename.single-c.freq", glob(catfile($single_c_dir,"*.freq")));

}

#######################################################################
# OLD countMethylation.pl

else{
    # Count methyl
    for my $base (sort values %base_gff_split) {
        my $singlec = catfile($single_c_dir, basename($base, ".gff")) . ".single-c.gff";

        if ($pm->start == 0){
            $logger->info("Processing $base");

            launch("perl -S countMethylation.pl --ref $opt_reference --gff $base --output ?? --freq $singlec.freq --sort -d $opt_di_nuc_freqs", 
                expected => $singlec, dryrun => $dry);

            my %split_by_context;
            if (!$dry){
                %split_by_context = GFF::Split::split_feature($singlec, @contexts);
                $logger->info("result of context split of $singlec: \n" . join "\n", values %split_by_context);
            }
            else {
                $logger->info("DRY: GFF::Split::split_feature($singlec, @contexts); ");
            }

            for my $singlec_context (values %split_by_context) {
                my $m = "$singlec_context.merged";
                launch("perl -S window_by_fixed.pl -o ?? $singlec_context", expected => $m, dryrun => $dry);

                unless ($opt_no_windowing){
                    # make window file name
                    my $windows_base = basename($singlec_context);
                    if ($windows_base !~ s/\.single-c-/.w$opt_window_size-/){
                        $logger->logdie("$singlec_context- naming screwed up?");
                    }
                    my $windows = catfile($windows_dir, $windows_base);

                    launch("perl -S window_by_fixed.pl -w $opt_window_size --reference $opt_reference --output ?? --no-skip $m", 
                        expected => $windows, dryrun => $dry);
                }
            }
            $pm->finish;
        }
    };

    $pm->wait_all_children;

    launch("perl -S collect-freqs.pl -o $basename.single-c.freq $single_c_dir", dryrun => $dry);

    # not sure if this works? 
    # chdir $single_c_dir;
    # for my $cont (@contexts) {
    #     my @files = glob("*$cont*");
    #     launch("single_c_concat.pl " . join(" ", @files), dryrun => $dry);
    # }
}

my @concatenated_single_c_files;
my $methyl_stats = "$basename.methyl-stats.txt";

for my $cont (sort @contexts) {
    $logger->info("concatenated $cont single-c files");
    my @files = glob(catfile($single_c_dir, "*$cont*"));
    push @concatenated_single_c_files, single_c_concat(@files);
}

if ($opt_stats){
    if (! -s $methyl_stats){
        launch("perl -S gff-methyl-stats.pl -o $methyl_stats -t @concatenated_single_c_files", dryrun => $dry);
    }
}

if ($opt_ends){
    my $threads = $opt_parallel * 2 || 1;
    launch("perl -S ends_analysis_batch.pl -a -p $threads -c $opt_ends -d $opt_out_directory", dryrun => $dry);
}

=head1 NAME

 bs-sequel.pl - Run dzlab's bisulfite sequencing analysis pipeline

=head1 SYNOPSIS

Simplest possible options for arabidopsis

 bs-sequel.pl -f /path/to/reference_genome.fasta -l /path/to/reads.fastq -b H2AZ_HO_6_BS -rnd 1

Single ends example. 

 bs-sequel.pl -f /path/to/genomes/TAIR_reference.fas -l reads.fastq -b bsseqrun -k 0 -n 2 -t Arabidopsis -ls 1 45 -rs 46 76 -2 0 -rnd 0 -mh 10 -d output dir

Paired ends example, where s_7_1_sequence.txt is the left read and s_7_2_sequence is the right.

 bs-sequel.pl -f /path/to/genomes/TAIR_reference.fas -l s_7_1_sequence.txt -r s_7_2_sequence.txt -b paired_ends -k 300 -n 2 -t leco --batch 1 -ls 1 40 -rs 41-51

Notes: The current recommended settings from Daniel are -mh 10 and -rnd 1. If
the right splice is long enough, such as 1-50 and 51-100 for a 100bp read, -2 1
is also recommend. Previously, for arabidopsis, we used '-rnd 1'.

=head1 DESCRIPTION

=head1 OPTIONS

=over 1

=item -f <fasta> | --reference <fasta> 

Reference genome file in Fasta format. Required.

=for Euclid
    fasta.type:        readable

=item -l <fastq> | --left-read <fastq>

Left reads file in fastq format. For single ends, only -l should be given. Required.

=for Euclid
    fastq.type:        readable

=item -r <fastq> | --right-read <fastq>

Left reads file in fastq format. For single ends, should not be given. Optional.

=for Euclid
    fastq.type:        readable

=item -ls <start> <end> | --left-splice <start> <end>

Start and end coordinate for the chunk of the --left-read fastq file to use for the left alignment.  
If you omit, full length will be used.
For example, if you are doing single ends with 1-45 and 46-76, use "-ls 1 45".

=item -rs <start> <end> | --right-splice <start> <end>

Start and end coordinate for the chunk of the --right-read fastq file to use for the right alignment.  Optional.
For example, if you are doing single ends with 1-45 and 46-76, use "-rs 46 76".

=item -b <label> | --base-name <label>

Label for the file names... can be anything.  Choose something descriptive. Required

=item -d <dir> | --out-directory <dir>

Directory to put all result files. If omitted, assumed to be the same as --base-name.

=item -k <len> | --library-size <len>

Approx length of the molecules in PAIRED ends. Default to 0 (which is the only valid value for single-ends).  
For paired ends, recommended is 300.  This parameter helps deal with possible inserts between the left and right reads:

  |------------------------------------| (original read, paired end)
  |----------->                          1', aligned to c2t
               |---------|               insert
                          <------------| 2', aligned to g2a

For single ends, there isn't an insert use 0 (default).

  |----------------------|               (original read, single end)
  |----------->                          1'
               |                         insert (size 0)
                |-------->               2' (simulated. compared against c2t just like 1')

=item -w <size> | --window-size <size>

Default 50, for windowing single-c files.

=for Euclid
    size.default:     50

=back

=head2 BOWTIE OPTIONS

=over 

=item -mh <hits> | --max-hits <hits>

Discards reads that map to the genome more the this many times, passed to bowtie.  In repetitive sections of the genome,
reads can potentially map hundreds of times, so this helps us filter repetitive chunks out..  Defaults to 0  for no filtering.  
Daniel says 10 is a good number to use.  Use 0 to disable. Default 10.

=for Euclid
    hits.default:     10

=item -n <num> | --mismatches <num>

For bowtie.  For each read alignment, allow this many mismatches.  Default 2.

=for Euclid
    num.default:     2

=item  -bf | --bowtie-file

Produce intermediate bowtie file.

=back

=head2 CORRELATION OPTIONS

=over 

=item -rnd <rnd> | --random-assign <rnd>

For correlatedPairedEnds.pl.  When there are multiple possible reconciliations between the left and right alignments,
assign one randomly.  For Arabidopsis and other organisms with lower levels of repetitive sequences, use 0. For maize
and similar, use 1.  Default 1.

=for Euclid
    rnd.default:     1

=item -2 <boolean> | --trust-dash-2 <boolean>

For correlatedPairedEnds.pl.  If 1, when the dash-2 reads align when the dash-1 reads don't, KEEP the read.  This
should be 0 for single ends since the downstream portion (right reads) of each read will be lower quality than the
upstream (left reads).  Default 0.

=for Euclid
    boolean.default:     0

=back

=head2 Less frequently used options

=over 

=item -t <orgname> | --organism <orgname>

Label for collect_align_stats.pl (the .log file produced).  Arabidopsis, Rice, Puffer, etc.

=for Euclid
    orgname.default:     'unknown'

=item -dnf <boolean> | --di-nuc-freqs <boolean> 

When this is 0, calculate the CG, CHH and CHG contexts.  If 1, calculate CG, CA, CT, CC.  Default 0.

=for Euclid
    boolean.default:     0

=item -i <num> | --batch <num>

Label for collect_align_stats.pl (the .log file produced).  Default to 1.

=for Euclid
    num.default:     1

=item --no-windowing

Skip the windowing after single-c file generation. 

=item --parallel <threads>

Number of subprocesses the script is allowed to run.  Default 0 (meaning no
subprocesses).  

=for Euclid
    threads.default:     0

=item --new-cm

Use the new countMethylation instead of the older countMethylation.pl. (No
longer experimental-- this is tested and ready to use).

=item  --dry 

=item  --downsample <fraction>

Downsample reads by given fraction.

=for Euclid
    fraction.type:        number, fraction >= 0 && fraction <= 1
    fraction.type.error:  <fraction> must be between 0 and 1.

=item  --ends <config>

Produces ends-analysis via ends_analysis_batch.pl with config file and --parallel * 2 threads.

=for Euclid
    config.type:        readable

=item --stats

Produces stats via gff-methyl-stats.pl

=item --verbose

Be verbose (only for --new-cm right now)

=item -h | --help

=back

=cut

