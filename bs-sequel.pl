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
use lib "$FindBin::Bin/lib";
use DZUtil qw/mfor timestamp split_names fastq_read_length/;
use Launch;
use GFF::Split;
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
# splice argument sanitizing

my $read_size = fastq_read_length($opt_left_read);
if (defined($opt_right_read) && $read_size =~ fastq_read_length($opt_right_read)){
    $logger->logdie("$opt_right_read and $opt_left_read not same read lengths?");
}

if (! %opt_left_splice){
    $logger->logwarn("Left splice -ls was not defined, using full length.");
    $opt_left_splice{start} = 1;
    $opt_left_splice{end} = $read_size;
}

my $do_right = %opt_right_splice;

if (! $do_right && !$opt_single_ends){
    $logger->logdie("You can't specify do paired ends without -rs!");
}

my @left_splice  = @opt_left_splice{qw/start end/};
my @right_splice = $do_right ? (@opt_right_splice{qw/start end/}) : ();

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
my $basename_right = catfile($opt_out_directory, basename($opt_left_read));
my $basename       = catfile($opt_out_directory, $opt_base_name);

#######################################################################
# Contexts

my @contexts;		# CG, CHG, CHH, etc.
unless (@contexts) {
    if ($opt_di_nuc_freqs) {@contexts = qw(CA CC CG CT)}
    else {@contexts = qw(CG CHG CHH)}
}

my $fasta_left = "$basename_left.fa";
my $fasta_right = "$basename_right.fa";

my $fasta_left_converted = "$basename_left.c2t";
my $fasta_right_converted = "$basename_right.g2a";

my $eland_left  = "${basename_left}_$left_splice[0]-$left_splice[1].eland3";
my $eland_right = $do_right ? "${basename_right}_$right_splice[0]-$right_splice[1].eland3" : "";

my $eland_left_post = "$eland_left.post";
my $eland_right_post = $do_right ? "$eland_right.post" : $eland_left_post;

if (! -f $eland_left_post || ($do_right && !-f $eland_right_post)){

#######################################################################
# convert reads

    launch("perl -S fq_all2std.pl fq2fa $opt_left_read > ??",  expected => $fasta_left, dryrun => $dry);
    launch("perl -S convert.pl c2t $fasta_left > ??", expected =>  $fasta_left_converted, dryrun => $dry);
    unless ($opt_single_ends) {
        launch("perl -S fq_all2std.pl fq2fa $opt_right_read > ??", expected =>  $fasta_right, dryrun => $dry);
        launch("perl -S convert.pl g2a $fasta_right > ??", expected =>  $fasta_right_converted, dryrun => $dry);
    }

#######################################################################
# convert genomes

    launch("perl -S rcfas.pl $opt_reference > ??", expected =>  "$opt_reference.rc", dryrun => $dry);
    launch("perl -S convert.pl c2t $opt_reference.rc > ??", expected =>  "$opt_reference.c2t", dryrun => $dry);
    unless ($opt_single_ends) {
        launch("perl -S convert.pl g2a $opt_reference.rc > ??", expected =>  "$opt_reference.g2a", dryrun => $dry);
    }

#######################################################################
# bowtie-build

    launch("bowtie-build $opt_reference.c2t $opt_reference.c2t", expected =>  "$opt_reference.c2t.1.ebwt", dryrun => $dry);
    unless ($opt_single_ends) {
        launch("bowtie-build $opt_reference.g2a $opt_reference.g2a", expected =>  "$opt_reference.g2a.1.ebwt", dryrun => $dry);
    }

#######################################################################
# bowtie

    my $l3trim = $read_size - $left_splice[1];
    my $l5trim = $left_splice[0] - 1;

    my $mh_args = $opt_max_hits ? " --strata  -k $opt_max_hits -m $opt_max_hits " : q{ };

# align with bowtie
    if ($pm->start == 0){
        launch("bowtie $opt_reference.c2t -f -B 1 -v $opt_mismatches -5 $l5trim -3 $l3trim --best $mh_args --norc $fasta_left_converted ??", 
            expected => $eland_left, dryrun => $dry, id => "left bowtie", accum => 1, also => $bowtie_logname);
        launch("perl -S parse_bowtie.pl -u $fasta_left -s @left_splice  $eland_left -o ??", 
            expected => $eland_left_post, dryrun => $dry, id => "left parse_bowtie");
        $pm->finish;
    }

    if ($pm->start == 0){
        if ($do_right){
            my $r3trim = $read_size - $right_splice[1];
            my $r5trim = $right_splice[0] - 1;
            if ($opt_single_ends) {
                launch("bowtie $opt_reference.c2t -f -B 1 -v $opt_mismatches -5 $r5trim -3 $r3trim --best $mh_args --norc $fasta_left_converted ??" , 
                    expected => $eland_right, dryrun => $dry, id => "right bowtie", accum => 1, also => $bowtie_logname);
                launch("perl -S parse_bowtie.pl -u $fasta_left -s @right_splice  $eland_right -o ??", 
                    expected => $eland_right_post, dryrun => $dry, id => "right parse_bowtie");
            }
            else {
                launch("bowtie $opt_reference.g2a -f -B 1 -v $opt_mismatches -5 $r5trim -3 $r3trim --best $mh_args --norc $fasta_right_converted ??" , 
                    expected => $eland_right, dryrun => $dry, id => "right bowtie", accum => 1, also => $bowtie_logname);
                launch("perl -S parse_bowtie.pl -u $fasta_right -s @right_splice  $eland_right -o ??", 
                    expected => $eland_right_post, dryrun => $dry, id => "right parse_bowtie");
            }
        }
        else {
            $eland_right_post = $eland_left_post;
        }
        $pm->finish;
    }

    $pm->wait_all_children;
}
else{
    $logger->info(".eland3.post already exists, don't re-create intermediate files up to then");
}

#######################################################################
# Correlate

my $base_gff = "$basename.gff";
my $base_log = "$basename.log";

if ($do_right){
    launch("perl -S correlatePairedEnds.pl -l $eland_left_post -r $eland_right_post -ref $opt_reference -o ?? -t 0 -d $opt_library_size -s $read_size -2 $opt_trust_dash_2 -1 $opt_single_ends -m $opt_max_hits -a $opt_random_assign", expected =>  $base_gff, dryrun => $dry);
}
else{
    launch("perl -S correlateSingleEnds.pl -e $eland_left_post --reference $opt_reference -o ?? -m $opt_max_hits --random-assign $opt_random_assign", expected =>  $base_gff, dryrun => $dry);
}

launch("perl -S collect_align_stats.pl $eland_left_post $eland_right_post $base_gff $opt_organism $opt_batch > ??", expected =>  $base_log, dryrun => $dry);

#######################################################################
# Split correlate

$logger->info("Splitting $base_gff by group");

my @base_gff_split=();
if (!$dry){
    @base_gff_split = GFF::Split::split_sequence($base_gff,@groups);
    $logger->info("result of context split of $base_gff: \n" . join "\n", @base_gff_split);
}
else {
    $logger->info("DRY: GFF::Split::split_sequence($base_gff,@groups);");
}

#######################################################################
# Count methyl

my @single_c_split = map {
    my $single_c_base = basename($_, ".gff");
    catfile($single_c_dir, $single_c_base) . ".single-c.gff";
} @base_gff_split;

mfor \@base_gff_split, \@single_c_split, sub{
    my ($base, $singlec) = @_;

    if ($pm->start == 0){
        $logger->info("Processing $base");

        launch("perl -S countMethylation.pl --ref $opt_reference --gff $base --output ?? --freq $singlec.freq --sort -d $opt_di_nuc_freqs", 
            expected => $singlec, dryrun => $dry, id => "count-$singlec");

        my @split_by_context = (); 
        if (!$dry){
            @split_by_context = GFF::Split::split_feature($singlec, @contexts);
            $logger->info("result of context split of $singlec: \n" . join "\n", @split_by_context);
        }
        else {
            $logger->info("DRY: GFF::Split::split_feature($singlec, @contexts); ");
        }

        for my $singlec_context (@split_by_context) {
            my $m = "$singlec_context.merged";
            #launch("perl -S compile_gff.pl -o ?? $singlec_context", expected => $m, dryrun => $dry, id => "compile-$singlec");
            launch("perl -S window_by_fixed.pl -o ?? $singlec_context", expected => $m, dryrun => $dry, id => "w1-$singlec");

            unless ($opt_no_windowing){
                # make window file name
                my $windows_base = basename($singlec_context);
                if ($windows_base !~ s/\.single-c-/.w$opt_window_size-/){
                    $logger->logdie("$singlec_context- naming screwed up?");
                }
                my $windows = catfile($windows_dir, $windows_base);

                launch("perl -S window_by_fixed.pl -w $opt_window_size --reference $opt_reference --output ?? --no-skip $m", 
                    expected => $windows, dryrun => $dry, id => "window-$singlec");
            }
        }
        $pm->finish;
    }
};


$pm->wait_all_children;

launch("perl -S collect-freqs.pl -o $basename.single-c.freq $single_c_dir", dryrun => $dry);

chdir $single_c_dir;
for my $cont (@contexts) {
    my @files = glob("*$cont*.merged");
    launch("single_c_concat.pl " . join(" ", @files), dryrun => $dry);
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

Launch various subprocesses (countMethylation only at the moment in parallel.)  Careful! 
Default 0 (meaning no subprocesses).

=for Euclid
    threads.default:     0

=item  --dry 

=back

=over

=item -h | --help

=back

=cut

