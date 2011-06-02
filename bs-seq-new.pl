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

use lib "$FindBin::Bin/lib";
use DZUtil qw/timestamp split_names/;
use Launch;
use GFF::Split;

my $logname = $opt_out_directory . "-" . timestamp() . ".log.txt";

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
    log4perl.appender.File.layout.ConversionPattern = %d{HH:mm:ss} %p.1> (%L) %m%n
/;
Log::Log4perl::init( \$conf );
my $logger = get_logger("PipeLine");

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
unless (
    $opt_left_read && $opt_right_read && $opt_reference && $opt_out_directory && $opt_base_name && $opt_read_size
);

my $pm = Parallel::ForkManager->new($opt_parallel);

my @left_splice = @opt_left_splice{qw/start end/};
my @right_splice = @opt_right_splice{qw/start end/};

my @contexts;		# CG, CHG, CHH, etc.
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

@left_splice  = (1, $opt_read_size) unless @left_splice;
@right_splice = (1, $opt_read_size) unless @right_splice;

#######################################################################
# Directory creation

my $windows_dir  = catfile ($opt_out_directory, 'windows');
my $single_c_dir = catfile ($opt_out_directory, 'single-c');

# check whether output directory should be created, exists, should overwritten
if (-d $opt_out_directory and $opt_overwrite) {rmtree ( $opt_out_directory, {keep_root => 1} )}
elsif (! -d $opt_out_directory) {
    mkpath ( $opt_out_directory ,  {verbose => 1} );
}
else {warn " overwrite $opt_out_directory"}

if (! -d $windows_dir) { mkpath ( $windows_dir,  {verbose => 1} ); }
if (! -d $single_c_dir){ mkpath ( $single_c_dir, {verbose => 1} ); }


unless (@contexts) {
    if ($opt_di_nuc_freqs) {@contexts = qw(CA CC CG CT)}
    else {@contexts = qw(CG CHG CHH)}
}

my $basename_left = catfile($opt_out_directory, basename($opt_left_read));
my $basename_right = catfile($opt_out_directory, basename($opt_left_read));

my %files = (
    # left fasta, and bisulfite treated
    lfa   => catfile ($opt_out_directory, basename($opt_left_read))  . '.fa',
    lc2t  => catfile ($opt_out_directory, basename($opt_left_read))  . '.c2t',
    # right fasta, and bisulfite treated
    rfa   => catfile ($opt_out_directory, basename($opt_right_read)) . '.fa',
    rg2a  => catfile ($opt_out_directory, basename($opt_right_read)) . '.g2a',
    # left/right eland3, from bowtie
    lel3  => catfile ($opt_out_directory, basename($opt_left_read))  . "_$left_splice[0]-$left_splice[1].eland3",
    rel3  => catfile ($opt_out_directory, basename($opt_right_read)) . "_$right_splice[0]-$right_splice[1].eland3",
    # output of correlated paired ends
    base  => catfile ($opt_out_directory, $opt_base_name)  . '.gff',
    # output of collect align states
    log   => catfile ($opt_out_directory, $opt_base_name)  . '.log',
    
    # <output dir>/basename-sequence.gff
    split => _gen_files (catfile ($opt_out_directory, $opt_base_name), 'gff',  @groups),

    # <output dir>/single-c/basename-sequence.single-c.gff
    freq  => _gen_files (catfile ($opt_out_directory, 'single-c', $opt_base_name), 'single-c.gff', @groups),

    # <output dir>/single-c/basename-sequence.single-c-CG.gff, 
    # <output dir>/single-c/basename-sequence.single-c-CHG.gff, etc.
    cont  => [map { _gen_files (catfile ($opt_out_directory, 'single-c', $opt_base_name), "single-c-$_.gff", @groups) } @contexts],

    # <output dir>/windows/basename-sequence.w50-CG.gff, 
    # <output dir>/windows/basename-sequence.w50-CHG.gff, etc.
    wcont => [map { _gen_files (catfile ($opt_out_directory, 'windows', $opt_base_name), "w${opt_window_size}-$_.gff", @groups) } @contexts],
);

#######################################################################
# convert reads

my $left_fasta = "$basename_left.fa";
my $right_fasta = "$basename_right.fa";

launch("perl -S fq_all2std.pl fq2fa $opt_left_read > ??",  expected => $files{lfa});
launch("perl -S convert.pl c2t $files{lfa} > ??", expected =>  $files{lc2t});
unless ($opt_single_ends) {
    launch("perl -S fq_all2std.pl fq2fa $opt_right_read > ??", expected =>  $files{rfa});
    launch("perl -S convert.pl g2a $files{rfa} > ??", expected =>  $files{rg2a});
}

#######################################################################
# convert genomes

launch("perl -S rcfas.pl $opt_reference > ??", expected =>  "$opt_reference.rc");
launch("perl -S convert.pl c2t $opt_reference.rc > ??", expected =>  "$opt_reference.c2t");
unless ($opt_single_ends) {
    launch("perl -S convert.pl g2a $opt_reference.rc > ??", expected =>  "$opt_reference.g2a");
}

#######################################################################
# align with bowtie

launch("bowtie-build $opt_reference.c2t $opt_reference.c2t", expected =>  "$opt_reference.c2t.1.ebwt");
unless ($opt_single_ends) {
    launch("bowtie-build $opt_reference.g2a $opt_reference.g2a", expected =>  "$opt_reference.g2a.1.ebwt");
}

my $l3trim = $opt_read_size - $left_splice[1];
my $l5trim = $left_splice[0] - 1;

my $r3trim = $opt_read_size - $right_splice[1];
my $r5trim = $right_splice[0] - 1;

my $mh_args = $opt_max_hits ? " --strata  -k $opt_max_hits -m $opt_max_hits " : q{ };
# align with bowtie
launch("bowtie $opt_reference.c2t -f -B 1 -v $opt_mismatches -5 $l5trim -3 $l3trim --best $mh_args --norc $files{lc2t} ??", expected => $files{lel3});
unless ($opt_single_ends) {
	launch("bowtie $opt_reference.g2a -f -B 1 -v $opt_mismatches -5 $r5trim -3 $r3trim --best $mh_args --norc $files{rg2a} ??" , expected => $files{rel3});
}
else {
	launch("bowtie $opt_reference.c2t -f -B 1 -v $opt_mismatches -5 $r5trim -3 $r3trim --best $mh_args --norc $files{lc2t} ??" , expected => $files{rel3});
}

#######################################################################
# parse bowtie

# get back original non-converted reads and convert from bowtie to eland3

launch("perl -S parse_bowtie.pl -u $files{lfa} -s @left_splice  $files{lel3} -o ??", expected => "$files{lel3}.post");
unless ($opt_single_ends) {
	launch("perl -S parse_bowtie.pl -u $files{rfa} -s @right_splice  $files{rel3} -o ??", expected => "$files{rel3}.post");
}
else {
	launch("perl -S parse_bowtie.pl -u $files{lfa} -s @right_splice  $files{rel3} -o ??", expected => "$files{rel3}.post");
}

# make sure reads map together
launch("perl -S correlatePairedEnds.pl -l $files{lel3}.post -r $files{rel3}.post -ref $opt_reference -o ?? -t 0 -d $opt_library_size -s $opt_read_size -2 $opt_trust_dash_2 -1 $opt_single_ends -m $opt_max_hits -a $opt_random_assign", expected =>  $files{base});

# basic stats about the aligment
launch("perl -S collect_align_stats.pl $files{lel3}.post $files{rel3}.post $files{base} $opt_organism $opt_batch > ??", expected =>  $files{log});

$logger->info("Splitting $files{base} by group");
GFF::Split::split_sequence($files{base},@groups);

# quantify methylation
for (@groups) {
    $pm->start and next;
    launch("perl -S countMethylation.pl --ref $opt_reference --gff $files{split}->{$_} --output $files{freq}->{$_} --sort -d $opt_di_nuc_freqs", expected => $files{freq}->{$_});
    $pm->finish;
}
$pm->wait_all_children;

# window methylation counts into non-overlapping windows
for my $context (0 .. @contexts - 1) {
    for my $group (@groups) {
        $logger->info("Splitting $files{freq}{$group} by context");
        GFF::Split::split_feature($files{freq}{$group}, @contexts);
        unless ($opt_no_windowing){
            my $m = "$files{cont}->[$context]{$group}.merged";
            launch("perl -S compile_gff.pl -o ?? $files{cont}->[$context]{$group}", expected => $m);
            launch("perl -S window_gff.pl $m --width $opt_window_size --step $opt_window_size --output ?? --no-skip",
            expected => $files{wcont}->[$context]{$group});
        }
    }
}

### DONE

sub _gen_files {
    my ($base, $ext, @groups) = @_;

    my %split_files;
    for my $group (@groups) {
        $split_files{$group} = "$base-$group.$ext";
    }
    return \%split_files;
}

=head1 NAME

 bs-seq-new.pl - Run dzlab's bisulfite sequencing analysis pipeline

=head1 SYNOPSIS

Paired ends example, where s_7_1_sequence.txt is the left read and s_7_2_sequence is the right.

  bs-seq-new.pl --left-read s_7_1_sequence.txt --right-read s_7_2_sequence.txt --reference REF_COL.fa --base-name test-out --read-size 45 --library-size 300 --mismatches 2 --organism leco --batch 1 --left-splice 1 40 --right-splice 5 45

Single ends example. Notice the left and right reads are the same file.

 bs-seq-new.pl -l reads.fastq -r reads.fastq -f /path/to/genomes/TAIR_reference.fas -b bsseqrun -s 76 -k 0 -n 2 -t Arabidopsis -i 1 -ls 1 45 -rs 46 76 -1 1 -2 0 -rnd 0 -mh 10 -d output dir

=head1 DESCRIPTION


=head1 OPTIONS

=over 1

=item -f <fasta> | --reference <fasta> 

Reference genome file in Fasta format. Required.

=for Euclid
    fasta.type:        readable

=item -l <fastq> | --left-read <fastq>

Left reads file in fastq format. For single ends, -l and -r should be the same. Required.

=for Euclid
    fastq.type:        readable

=item -r <fastq> | --right-read <fastq>

Left reads file in fastq format. For single ends, -l and -r should be the same. Required.

=for Euclid
    fastq.type:        readable

=item -ls <start> <end> | --left-splice <start> <end>

Start and end coordinate for the chunk of the --left-read fastq file to use for the left alignment.  Required.
For example, if you are doing single ends with 1-45 and 46-76, use "-ls 1 45".

=item -rs <start> <end> | --right-splice <start> <end>

Start and end coordinate for the chunk of the --right-read fastq file to use for the right alignment.  Required.
For example, if you are doing single ends with 1-45 and 46-76, use "-rs 46 76".

=item -s <len> | --read-size <len> 

Length of each read.  76, 100, etc. Required.

=for Euclid
    len.default:     76

=item -b <label> | --base-name <label>

Label for the file names... can be anything.  Choose something descriptive. Required

=item -d <dir> | --out-directory <dir>

Directory to put all result files. Required.

=item -o | --overwrite 

Enable to overwrite the --directory if it already exists.  

=item -k <len> | --library-size <len>

Approx length of the molecules in PAIRED ends. Default 300.  For single ends, use 0.  Ie, for paired ends, this
parameter helps deal with possible inserts between the left and right reads:

  |------------------------------------| (original read, paired end)
  |----------->                          1', aligned to c2t
               |---------|               insert
                          <------------| 2', aligned to g2a

For single ends, there isn't an insert use 0.

  |----------------------|               (original read, single end)
  |----------->                          1'
               |                         insert (size 0)
                |-------->               2' (simulated. compared against c2t just like 1')

=for Euclid
    len.default:     300

=item -t <orgname> | --organism <orgname>

Label for collect_align_stats.pl (the .log file produced).  Arabidopsis, Rice, Puffer, etc.

=for Euclid
    orgname.default:     '.'

=item -w <size> | --window-size <size>

Default 50, for windowing single-c files.

=for Euclid
    size.default:     50

=item -1 <boolean> | --single-ends <boolean>

1 if single ends, 0 if paired.  Default 1.

=for Euclid
    boolean.default:     1

=back

=head2 BOWTIE OPTIONS

=over 

=item -mh <hits> | --max-hits <hits>

Discards reads that map to the genome more the this many times, passed to bowtie.  In repetitive sections of the genome,
reads can potentially map hundreds of times, so this helps us filter repetitive chunks out..  Defaults to 0  for no filtering.  
Daniel says 10 is a good number to use.  Use 0 (default) to disable.

=for Euclid
    hits.default:     0

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

=back

=over

=item -h | --help

=back

=cut

