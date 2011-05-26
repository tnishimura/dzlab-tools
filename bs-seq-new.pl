#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use File::Spec;
use File::Path;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/lib";
use DZUtil qw/split_names/;


# SIGINT trap. Ctrl-c triggers unclean exit.
$SIG{INT} = sub {croak "Received SIG$_[0]. Exiting...\n"};

my $left_read;   # required
my $right_read;  # required
my $reference;   # required
my $base_name     = 'out';
my $overwrite;
my $read_size     = 76;     # length of the reads
my $library_size  = 300;    # for paired ends, approx size of the reads.  For single, use 0
my $mismatches    = 2;      # num mismatches to allow for bowtie
my $organism      = q{.};   # label for collectStats.
my $batch         = 1;      # label for collectStats.
my $window_size   = 50;     
my $trust_dash_2  = 0;
my $single_ends   = 1;
my @left_splice;
my @right_splice;
my @groups        = (); 	# sequences (chr1, chr2, ...)
my $aligner       = 'bowtie';
my $max_hits      = 0;      # cutoff for bowtie alighnment before discarding
my $random_assign = 1;      # if there is multiple reconciliation between ends, use one randomly.  Otherwise discard.
my $pthreads      = 1;      # not used... for bowtie.
my $di_nuc_freqs  = 0;      # if true, CG, CT, CC, CA instead of CHH, CHG, CG. 
my @contexts;		# CG, CHG, CHH, etc.
my $no_windowing  = 0;
my $pfm_max       = 1;

my @date = localtime (time);

my $out_dir = File::Spec->catdir (
    File::Spec->curdir(),
    sprintf "DZ_full-run_%4d-%02d-%02d_%02d.%02d.%02d",
    $date[5] + 1900, $date[4] + 1, $date[3], $date[2], $date[1], $date[0]
);

# Grabs and parses command line options
my $result = GetOptions (
    'left-read|l=s'        => \$left_read,
    'right-read|r=s'       => \$right_read,
    'reference|f=s'        => \$reference,
    'base-name|b=s'        => \$base_name,
    'overwrite|o'          => \$overwrite,
    'read-size|s=i'        => \$read_size,
    'library-size|k=i'     => \$library_size,
    'mismatches|n=i'       => \$mismatches,
    'organism|t=s'         => \$organism,
    'batch|i=i'            => \$batch,
    'window-size|w=i'      => \$window_size,
    'trust-dash-2|2=i'     => \$trust_dash_2,
    'single-ends|1=i'      => \$single_ends,
    'left-splice|ls=i{2}'  => \@left_splice,
    'right-splice|rs=i{2}' => \@right_splice,
    'groups|g=s{,}'        => \@groups,
    'out-directory|d=s'    => \$out_dir,
    'aligner|a=s'          => \$aligner,
    'max-hits|mh=i'        => \$max_hits,
    'random-assign|rnd=i'  => \$random_assign,
    'pthreads|pt=i'        => \$pthreads,
    'di-nuc-freqs|dnf=i'   => \$di_nuc_freqs,
    'contexts|ct=s{,}'     => \@contexts,
    'no-windowing'         => \$no_windowing,
    'parallel=i'           => \$pfm_max,
    'verbose|V'            => sub { use diagnostics; },
    'quiet|q'              => sub { no warnings; },
    'help|h'               => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'             => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless $result and $left_read and $right_read and $reference;

use Parallel::ForkManager;
my $pm = Parallel::ForkManager->new($pfm_max);

# Get all chromosomes, pseudo-chromosomes, groups, etc, in the fasta reference file
# Discards all information after first blank character in fasta header
unless (@groups) {
    open my $REFERENCE, '<', $reference or croak "Can't open $reference: $!";

    while (<$REFERENCE>) {
        if (m/>/) {
            tr/A-Z/a-z/;
            m/>([^\s]+)/ && push @groups, $1;
        }
    }
    close $REFERENCE or carp "Can't close $reference: $!";
}

@left_splice  = (1, $read_size) unless @left_splice;
@right_splice = (1, $read_size) unless @right_splice;

# check whether output directory should be created, exists, should overwritten
if (-d $out_dir and $overwrite) {rmtree ( $out_dir, {keep_root => 1} )}
elsif (! -d $out_dir) {
    mkpath ( File::Spec->catfile ($out_dir, 'windows'),  {verbose => 1} );
    mkpath ( File::Spec->catfile ($out_dir, 'single-c'), {verbose => 1} );
}
else {warn " overwrite $out_dir"}

unless (@contexts) {
    if ($di_nuc_freqs) {@contexts = qw(CA CC CG CT)}
    else {@contexts = qw(CG CHG CHH)}
}

my %files = (
    # left fasta, and bisulfite treated
    lfa   => File::Spec->catfile ($out_dir, basename($left_read))  . '.fa',
    lc2t  => File::Spec->catfile ($out_dir, basename($left_read))  . '.c2t',
    # right fasta, and bisulfite treated
    rfa   => File::Spec->catfile ($out_dir, basename($right_read)) . '.fa',
    rg2a  => File::Spec->catfile ($out_dir, basename($right_read)) . '.g2a',
    # left/right eland3, from bowtie
    lel3  => File::Spec->catfile ($out_dir, basename($left_read))  . "_$left_splice[0]-$left_splice[1].eland3",
    rel3  => File::Spec->catfile ($out_dir, basename($right_read)) . "_$right_splice[0]-$right_splice[1].eland3",
    # output of correlated paired ends
    base  => File::Spec->catfile ($out_dir, $base_name)  . '.gff',
    # output of collect align states
    log   => File::Spec->catfile ($out_dir, $base_name)  . '.log',
    
    # <output dir>/basename-sequence.gff
    split => _gen_files (File::Spec->catfile ($out_dir, $base_name), 'gff',  @groups),

    # <output dir>/single-c/basename-sequence.single-c.gff
    freq  => _gen_files (File::Spec->catfile ($out_dir, 'single-c', $base_name), 'single-c.gff', @groups),

    # <output dir>/single-c/basename-sequence.single-c-CG.gff, 
    # <output dir>/single-c/basename-sequence.single-c-CHG.gff, etc.
    cont  => [map { _gen_files (File::Spec->catfile ($out_dir, 'single-c', $base_name), "single-c-$_.gff", @groups) } @contexts],

    # <output dir>/windows/basename-sequence.w50-CG.gff, 
    # <output dir>/windows/basename-sequence.w50-CHG.gff, etc.
    wcont => [map { _gen_files (File::Spec->catfile ($out_dir, 'windows', $base_name), "w${window_size}-$_.gff", @groups) } @contexts],
);

# convert reads
run_cmd ("perl -S fq_all2std.pl fq2fa $left_read > $files{lfa}")  unless file_exists($files{lfa});
run_cmd ("perl -S convert.pl c2t $files{lfa} > $files{lc2t}")     unless file_exists($files{lc2t});
unless ($single_ends) {
    run_cmd ("perl -S fq_all2std.pl fq2fa $right_read > $files{rfa}") unless file_exists($files{rfa});
    run_cmd ("perl -S convert.pl g2a $files{rfa} > $files{rg2a}")     unless file_exists($files{rg2a});
}

# convert genomes
run_cmd ("perl -S rcfas.pl $reference > $reference.rc")           unless file_exists("$reference.rc");
run_cmd ("perl -S convert.pl c2t $reference.rc > $reference.c2t") unless file_exists("$reference.c2t");
run_cmd ("perl -S convert.pl g2a $reference.rc > $reference.g2a") unless file_exists("$reference.g2a") or $single_ends;

# align with bowtie
run_cmd ("bowtie-build $reference.c2t $reference.c2t") unless file_exists("$reference.c2t.1.ebwt");
run_cmd ("bowtie-build $reference.g2a $reference.g2a") unless file_exists("$reference.g2a.1.ebwt") or $single_ends;

my $l3trim = $read_size - $left_splice[1];
my $l5trim = $left_splice[0] - 1;

my $r3trim = $read_size - $right_splice[1];
my $r5trim = $right_splice[0] - 1;

# align with bowtie
run_cmd ("bowtie $reference.c2t -f -B 1 -v $mismatches -5 $l5trim -3 $l3trim --best" . ($max_hits ? " --strata  -k $max_hits -m $max_hits" : q{}) . " --norc $files{lc2t} $files{lel3}" )   unless file_exists($files{lel3});
unless ($single_ends) {
	run_cmd ("bowtie $reference.g2a -f -B 1 -v $mismatches -5 $r5trim -3 $r3trim --best" . ($max_hits ? "  --strata -k $max_hits -m $max_hits" : q{}) . " --norc $files{rg2a} $files{rel3}" ) unless file_exists($files{rel3});
}
else {
	run_cmd ("bowtie $reference.c2t -f -B 1 -v $mismatches -5 $r5trim -3 $r3trim --best" . ($max_hits ? "  --strata -k $max_hits -m $max_hits" : q{}) . " --norc $files{lc2t} $files{rel3}" ) unless file_exists($files{rel3});
}

# get back original non-converted reads and convert from bowtie to eland3
run_cmd ("perl -S parse_bowtie.pl -u $files{lfa} -s @left_splice  $files{lel3} -o $files{lel3}.post") unless file_exists("$files{lel3}.post");
unless ($single_ends) {
	run_cmd ("perl -S parse_bowtie.pl -u $files{rfa} -s @right_splice  $files{rel3} -o $files{rel3}.post") unless file_exists("$files{rel3}.post");
}
else {
	run_cmd ("perl -S parse_bowtie.pl -u $files{lfa} -s @right_splice  $files{rel3} -o $files{rel3}.post") unless file_exists("$files{rel3}.post");
}

# make sure reads map together
run_cmd ("perl -S correlatePairedEnds.pl -l $files{lel3}.post -r $files{rel3}.post -ref $reference -o $files{base} -t 0 -d $library_size -s $read_size -2 $trust_dash_2 -1 $single_ends -m $max_hits -a $random_assign") unless file_exists($files{base});

# basic stats about the aligment
run_cmd ("perl -S collect_align_stats.pl $files{lel3}.post $files{rel3}.post $files{base} $organism $batch > $files{log}") unless file_exists($files{log});



# quantify methylation
for (@groups) {
    run_cmd ("perl -S split_gff.pl --sequence all $files{base}") unless (file_exists($files{split}->{$_}));
    $pm->start and next;
    run_cmd ("perl -S countMethylation.pl --ref $reference --gff $files{split}->{$_} --output $files{freq}->{$_} --sort -d $di_nuc_freqs") unless file_exists($files{freq}->{$_});
    $pm->finish;
}
$pm->wait_all_children;



# window methylation counts into non-overlapping windows
for my $context (0 .. @contexts - 1) {
    for my $group (@groups) {
        run_cmd ("perl -S split_gff.pl --feature all $files{freq}->{$group}") 
		unless file_exists($files{cont}->[$context]{$group});
        unless ($no_windowing){
            run_cmd ("perl -S window_gff.pl $files{cont}->[$context]{$group} --width 1 --step 1 --output $files{cont}->[$context]{$group}.merged") 
            unless file_exists("$files{cont}->[$context]{$group}.merged");
            run_cmd ("perl -S window_gff.pl $files{cont}->[$context]{$group}.merged --width $window_size --step $window_size --output $files{wcont}->[$context]{$group} --no-skip") 
            unless file_exists($files{wcont}->[$context]{$group});
        }
    }
}

### DONE


sub _error {
    my ($e) = @_;

    if ($e == -1) {
        return "failed to execute: $!\n";
    }
    elsif ($e & 127) {
        return sprintf "child died with signal %d, %s coredump\n",
        ($e & 127),  ($e & 128) ? 'with' : 'without';
    }
    else {
        return sprintf "child exited with value %d\n", $e >> 8;
    }
}

sub run_cmd {
    my ($cmd) = @_;
    warn "-- CMD: $cmd\n";
    eval {system ("$cmd") == 0 or croak _error( $? )};
    croak "** failed to run command '$cmd': $@" if $@;
}

sub _gen_files {
    my ($base, $ext, @groups) = @_;

    my %split_files;
    for my $group (@groups) {
        $split_files{$group} = "$base-$group.$ext";
    }
    return \%split_files;
}

sub file_exists {
    my $file = shift;
    print "$file exists: skipping..." && return 1 if -f $file and -s $file;
    return 0;
}



=head1 NAME

 bs-seq.pl - Run dzlab's bisulfite sequencing analysis pipeline

=head1 SYNOPSIS

Paired ends example, where s_7_1_sequence.txt is the left read and s_7_2_sequence is the right.

  bs-seq.pl --left-read s_7_1_sequence.txt --right-read s_7_2_sequence.txt --reference REF_COL.fa --base-name test-out --read-size 45 --library-size 300 --mismatches 2 --organism leco --batch 1 --left-splice 1 40 --right-splice 5 45

Single ends example. Notice the left and right reads are the same file.

 bs-seq.pl -l reads.fastq -r reads.fastq -f /path/to/genomes/TAIR_reference.fas -b bsseqrun -s 76 -k 0 -n 2 -t Arabidopsis -i 1 -ls 1 45 -rs 46 76 -1 1 -2 0 -rnd 0 -mh 10 -d output dir

=head1 DESCRIPTION


=head1 OPTIONS

=over 1

=item -f <fasta> | --reference <fasta> 

Reference genome file in Fasta format. Required.

=item -l <fastq> | --left-read <fastq>

Left reads file in fastq format. For single ends, -l and -r should be the same. Required.

=item -r <fastq> | --right-read <fastq>

Left reads file in fastq format. For single ends, -l and -r should be the same. Required.

=item -ls <start> <end> | --left-splice <start> <end>

Start and end coordinate for the chunk of the --left-read fastq file to use for the left alignment.  Required.
For example, if you are doing single ends with 1-45 and 46-76, use "-ls 1 45".

=item -rs <start> <end> | --right-splice <start> <end>

Start and end coordinate for the chunk of the --right-read fastq file to use for the right alignment.  Required.
For example, if you are doing single ends with 1-45 and 46-76, use "-rs 46 76".

=item -s <len> | --read-size <len> 

Length of each read.  76, 100, etc. Required.

=item -b <label> | --base-name <label>

Label for the file names... can be anything.  Choose something descriptive. Required

=item -d | --out-directory

Directory to put all result files. Required.

=item -o <boolean> | --overwrite <boolean>

Whether to overwrite the --directory if it already exists.  Default 0.

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

=item -t <orgname> | --organism <orgname>

Label for collect_align_stats.pl (the .log file produced).  Arabidopsis, Rice, Puffer, etc.

=item -w <size> | --window-size <size>

Default 50, for windowing single-c files.

=item -1 <boolean> | --single-ends <boolean>

1 if single ends, 0 if paired.  Default 1.

=back

=head2 BOWTIE OPTIONS

=over 

=item -mh <hits> | --max-hits <hits>

Discards reads that map to the genome more the this many times, passed to bowtie.  In repetitive sections of the genome,
reads can potentially map hundreds of times, so this helps us filter repetitive chunks out..  Defaults to 0  for no filtering.  
Daniel says 10 is a good number to use. 

=item -n <num> | --mismatches <num>

For bowtie.  For each read alignment, allow this many mismatches.  Default 2.

=back

=head2 Ends CORRELATION OPTIONS

=over 

=item -rnd | --random-assign

For correlatedPairedEnds.pl.  When there are multiple possible reconciliations between the left and right alignments,
assign one randomly.  For Arabidopsis and other organisms with lower levels of repetitive sequences, use 0. For maize
and similar, use 1.  Default 1.

=item -2 <boolean> | --trust-dash-2 <boolean>

For correlatedPairedEnds.pl.  If 1, when the dash-2 reads align when the dash-1 reads don't, KEEP the read.  This
should be 0 for single ends since the downstream portion (right reads) of each read will be lower quality than the
upstream (left reads).  Default 0.

=back

=head2 Less frequently used options

=over 

=item -dnf | --di-nuc-freqs

When this is 0, calculate the CG, CHH and CHG contexts.  If 1, calculate CG, CA, CT, CC.  Default 0.

=item -i <num> | --batch <num>

Label for collect_align_stats.pl (the .log file produced).  Default to 1.

=item --no-windowing

Skip the windowing after single-c file generation. 

=item --parallel <threads>

Launch various subprocesses (countMethylation only at the moment in parallel.)  Careful! 

=back

=over

=item -V | --verbose

=item -q | --quiet

=item -h | --help

=item -m | --manual

=back

=head1 REVISION

 Version 0.0.1

 $Rev: 440 $:
 $Author: psilva $:
 $Date: 2010-11-11 16:15:45 -0800 (Thu, 11 Nov 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/bs-seq.pl $:
 $Id: bs-seq.pl 440 2010-11-12 00:15:45Z psilva $:

=head1 AUTHOR

 Pedro Silva <pedros@berkeley.edu/>
 Zilberman Lab <http://dzlab.pmb.berkeley.edu/>
 Plant and Microbial Biology Department
 College of Natural Resources
 University of California, Berkeley

=head1 COPYRIGHT

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

