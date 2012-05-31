#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use File::Basename qw/basename/;
use File::Path qw/make_path remove_tree/;
use File::Spec::Functions qw/rel2abs catdir catfile/;
use Getopt::Long;
use List::MoreUtils qw/ any notall/;
use List::Util qw/sum/;
use Parallel::ForkManager;
use Pod::Usage;
use YAML qw/LoadFile DumpFile/;

use FindBin;
use lib "$ENV{HOME}/dzlab-tools/lib";
use FastqReader;
use FastaReader;
use DZUtil qw/c2t/;
use Digest::MD5::Util;
use Launch qw/cast/;
use Run::BowtieBuild;

END {close STDOUT}
$| = 1;

my @flank_range;

my $result = GetOptions (
    #"tdna-file|tf=s"       => \(my $tdna_file),
    #"reads-file|rf=s"      => \(my $reads_file),
    #"flank-seq-file|ff=s"  => \(my $flank_seq_file),
    "outdir|o=s"             => \(my $outdir),               # should already exist
    "reference|r=s"          => \(my $reference),            # reference file

    "flank-range|fr=i{2}"    => \@flank_range,               # this range of flank becomes range
    "flank-mismatch|fm=i"    => \(my $flank_mismatch = 3),   # align bait to reads with this many mm

    "prefix-trim5|p5=i"      => \(my $prefix_trim5 = 0),     # trim this much from prefix 5'
    "prefix-trim3|p3=i"      => \(my $prefix_trim3 = 5),     # trim this much from prefix 3'
    "prefix-min-size|pmin=i" => \(my $prefix_min_size = 20), # if trimmed is smaller, discard
    "prefix-max-size|pmax=i" => \(my $prefix_max_size = 30), # if trimmed is bigger, trim 3' further
    "prefix-mismatch|tm=i"   => \(my $prefix_mismatch = 3),  # align prefixes to tdna/ref with this many mm

    "parallel|p=i"           => \(my $parallel = 0),
);


@flank_range = @flank_range ? @flank_range : (1, 30);

usage("malformed arguments") if (!$result);
usage() if ! @ARGV;

my ($tdna_file, $flank_seq_file, $reads_file) = @ARGV;
usage("tdna, flank, reads file all need to exist") if any { ! defined $_ || ! -f $_ } ($tdna_file, $flank_seq_file, $reads_file);
usage("no reference given") if (!$reference);
usage("output directory dne") if (! defined $outdir || ! -d $outdir);

my $pm = Parallel::ForkManager->new($parallel);

sub usage { pod2usage(scalar(@_) ? (-msg => shift()) : (), -verbose => 2, -noperldoc => 1); }

#######################################################################

my $bait_file = catfile($outdir, "bait_$flank_range[0]_$flank_range[1]");
my $split_dir = catdir($outdir, "split_reads");
make_path $split_dir;
my $split_prefix = catdir($split_dir, basename($reads_file)) . ".";
my $split_checksum_file = $split_prefix . "MD5SUM";

my $tdna_file_bsrc; 
my $reference_bsrc;

#######################################################################
# create bait file from flank seq 
{
    my $fr = FastaReader->new(file => $flank_seq_file);
    my @seqs_in_flank = $fr->sequence_list();
    if (! @seqs_in_flank) {
        die "flank file empty";
    }
    open my $fh, '>', $bait_file;
    for my $flankseq (@seqs_in_flank) {
        my $bait = $fr->get($flankseq, @flank_range, bs => 'c2t', lenient => 1);
        say $fh ">${flankseq}_$flank_range[0]_$flank_range[1]";
        say $fh $bait;
    }
    close $fh;
}

#######################################################################
# create bsrc tdna

($tdna_file_bsrc) = bowtie_build(file => $tdna_file, rc => 1, bs => 'c2t');
($reference_bsrc) = bowtie_build(file => $reference, rc => 1, bs => 'c2t');

#######################################################################
# split the reads in smaller chunks (while converting to fasta/c2t)

my @split_read_files;
SPLIT:
{
    if (-f $split_checksum_file && defined(my $confirmed_files =md5_confirm($split_checksum_file))){
        @split_read_files = @$confirmed_files;
        say "splitting of reads was already done";
        last SPLIT;
    }

    my $max_per_file = 943_718_400; # 900MB
    my $suffix = 'aa';
    my $bytes_written = 0;;

    my $fqr = FastqReader->new(file => $reads_file, linesper => 4);
    open my $writer, '>', $split_prefix . $suffix;
    push @split_read_files, $split_prefix . $suffix;
    while (defined(my $fq = $fqr->next())){
        my ($readid, $sequence) = @$fq;
        $sequence = c2t $sequence;

        my $to_write = ">$readid\n$sequence\n";
        print $writer $to_write;

        $bytes_written += length $to_write;

        if ($bytes_written > $max_per_file){
            close $writer;
            $suffix++;
            open $writer, '>', $split_prefix . $suffix;
            push @split_read_files, $split_prefix . $suffix;
            $bytes_written = 0;
        }
    }
    close $writer;
    
    @split_read_files = map { rel2abs($_) } @split_read_files;

    md5_create_checksum($split_checksum_file, @split_read_files);
}

for my $split (@split_read_files) {
    $pm->start and next;
    bowtie_build(file => $split, noref => 1, rc => 0);
    $pm->finish; 
}
$pm->wait_all_children;

#######################################################################

my %aligned_reads; # split_read_file => { read_id => position }

# run each bowtie in separate process, collect reads and start positions
$pm->run_on_finish(sub{ 
        my $ref = $_[5];
        if (defined($ref)) {  
            my ($split_read_file, $id_to_startpos_href) = @$ref;
            $aligned_reads{$split_read_file} = $id_to_startpos_href;
        } 
    });

# each process returns [split_read_file, { readid => pos_in_readid, ... }]
for my $split_read (@split_read_files){
    my $bowtie_output = catfile($outdir, basename($bait_file) . "-vs-" . basename($split_read) . ".mm${flank_mismatch}.bowtie");

    $pm->start and next;
    if (! -f $bowtie_output){
        cast "bowtie --norc -f -v $flank_mismatch $split_read $bait_file > $bowtie_output",
    }
    else{
        say STDERR "bowtie already done ($bowtie_output)";
    }

    my $return_value = [ $split_read, {} ];

    # read the 
    open my $fh, '<', $bowtie_output;
    while (defined(my $line = <$fh>)){
        chomp $line;
        my (undef, undef, $readid, $position) = split /\t/, $line;
        $return_value->[1]{$readid} = $position;
    }
    close $fh;

    $pm->finish(0, $return_value);
}
$pm->wait_all_children;

#######################################################################

my $prefix_file = 
catfile($outdir, basename($bait_file)) . 
"-vs-" . 
basename($reads_file) . 
".mm${flank_mismatch}.prefix-trim$prefix_trim5-$prefix_trim3.size$prefix_min_size-$prefix_max_size.fasta";

if (! -f $prefix_file){
    say STDERR "grabbing prefixes from bowtie results";

    open my $prefix_out, '>', $prefix_file;
    while (my ($split_reads_file,$id_to_startpos) = each %aligned_reads) {
        next if (keys %$id_to_startpos == 0);

        # get original reads, since bowtie only keeps matching portion.
        my $fqr = FastqReader->new(file => $split_reads_file, linesper => 2);
        my $id2seq = $fqr->get_reads(keys %$id_to_startpos);


        while (my ($id,$seq) = each %$id2seq) {
            my $whole_prefix_size = $id_to_startpos->{$id}; # position reported in bait file bowtie is size of prefix
            my $length = $whole_prefix_size - $prefix_trim3 - $prefix_trim5;

            if ($length >= $prefix_min_size){
                my $offset = $prefix_trim5;
                
                if ($length > $prefix_max_size){
                    $length = $prefix_max_size;
                }

                {
                    my $start = $offset + 1;
                    my $end = $start + $length - 1;
                    say $prefix_out ">PREFIX_${id}_${start}_${end}";
                }

                say $prefix_out substr $seq, $offset, $length;;
            }
        }
    }
    close $prefix_out;
}
else{
    say STDERR "prefix file already done";
}

my $prefix_vs_tdna = "$prefix_file-vs-tdna.mm$prefix_mismatch";
my $prefix_vs_reference = "$prefix_file-vs-reference.mm$prefix_mismatch";

cast("bowtie --norc -f -B 1 -v $prefix_mismatch $tdna_file_bsrc $prefix_file > $prefix_vs_tdna");
cast("bowtie --norc -f -B 1 -v $prefix_mismatch $reference_bsrc $prefix_file > $prefix_vs_reference");

=head1 USAGE

 red-headed-step-child.pl [options] -r reference.fasta tdna.fasta flanks.fasta reads.fasta 

=head1 OPTIONS

=over

=item --outdir <dir> | -o <dir>          
 
Output directory. Required.

=item --reference <fasta>         | -r <fasta>        
 
Reference genome fasta file. Required.

=item --flank-range <start> <end> | -fr <start> <end> 

This range of flank becomes the bait (see below). Default 1, 30.

=item --flank-mismatch <n>        | -fm <n> 

Number of mismatches to use when bowtie-ing baits against reads. (default 3)

=item --prefix-trim5 <n>          | -p5 <n>

Number of bases to trim from 5' of the prefix of bait alignment in reads. Default 0.

=item --prefix-trim3 <n>          | -p3 <n> 

Number of bases to trim from 3' of the prefix of bait alignment in reads. Default 5.

=item --prefix-min-size <n>       | -pmin <n>      

If prefix, after trimming (as specified by options above) are smaller than this, discard.  Default 20.

=item --prefix-max-size <n>       | -pmax <n>       

If prefix, after trimming (as specified by options above) are larger than this, trim 3' further.  Default 30.

=item --prefix-mismatch <n>       | -tm <n>           

Number of mismatches to use when bowtie-ing prefixes against genome/tdna. (default 3)

=item --parallel <n>              | -p <n>

Run with this many threads.

=back

=head2 What it does

This script requires 4 files to run: the reference.fasta genome of the organism, tdna.fasta
containing the sequence of the tdna insertion, flanks.fasta containing the downstream, and reads.fasta 
containing the bisulfite-seq short reads. Schematically:

                        |----> tdna
                         \  /
                          \/
 |---------------------------------------------------------| genome
                           |---------> flank
                    |------------| reads that we want to find

=head2 How it works

=over

=item 1. Create baits from flank

First, we create baits from the flanks.  Baits are simply smaller subsections
of the flanking sequences.  One bait is created from each flank by taking the
subsection specified by --flank-range (default 1 to 30, meaning the baits will
be the first 30 bp of the flanks).  If the flank is not long enough, as much as
possible will be used. (If you're doing --flank-range 5 30 and the flank is
only 25bp, 5-25 of the flank will be the bait).

         bait == --flank-range
       <----------------------->
 5' |---------------------------------------------------> flank

=item 2. Build bowtie indices of the reads.

We are going to need to align the baits against the reads, and for that we need
to create bowtie indices (*.ebwt files) with bowtie-build.  

This is the opposite of the usual usage of bowtie. 

Note that the script splits reads into smaller (900 MB) chunks b/c bowtie
doesn't handle anything larger than 4gb. 

=item 3. Align baits against reads

Baits are aligned against reads with --flank-mismatch mismatches (default 3)
and --norc (so that only forward strand of reads are considered).  The reads found are
the ones which potentially mapped to junction of tdna insertion.

=item 4. Grab prefix of baited reads.

For each read baited, grab the prefix (the bases BEFORE where the bait aligns).
Mrim -p5 and -p3 from prefix (default 0 and 5).  If the prefix smaller than
-pmin (default 20), discard. If bigger than -pmax, further trim the 3' end. 

Schematically: 

  -p5             -p3
  <-->            <--> 
 |--------------------|------------------|---------------| Baited Read
      |-----------|           BAIT
          Prefix              

=item 5. Align prefixes against the tdna, references

The tdna and reference genome fasta file are both c2t converted.  Then, prefixes
are aligned against each with --prefix-mismatch mismatches (default 3). 

=item 6. You open the bowtie files from step 5, and count them.

Many hits in the -vs-tdna bowtie is evidence that the tdna was inserted.  Many
hits in -vs-reference in addition to -vs-tdna is evidence of heterozygous
insertion.

=back


=cut
