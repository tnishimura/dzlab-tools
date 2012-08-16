#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use File::Basename qw/basename/;
use File::Path qw/make_path/;
use File::Spec::Functions qw/catdir catfile/;
use Getopt::Long qw/:config no_ignore_case/;
use List::MoreUtils qw/ any notall/;
use List::Util qw/sum/;
use Pod::Usage;

use FindBin;
use lib "$ENV{HOME}/dzlab-tools/lib";

use FastqReader;
use DZUtil qw/c2t fastq_read_length/;
use Digest::MD5::Util;
use Launch qw/cast/;
use Run::BowtieBuild;
use Run::Bowtie2;

END {close STDOUT}
$| = 1;

my @splice;

my $result = GetOptions (
    "output-directory|d=s" => \(my $output_directory),
    "flank-prefix|fp=i"    => \(my $flank_prefix_length),
    "tdna-prefix|tp=i"     => \(my $tdna_prefix_length),

    "splice|s=i{2}"        => \@splice,

    # 4 necessary files
    "reference-file|r=s"   => \(my $reference_file),            # reference file
    "tdna-file|T=s"        => \(my $tdna_file),
    "flank-file|F=s"       => \(my $flank_file),
    "reads-file|R=s"       => \(my $reads_file),
);

usage("malformed arguments") if (!$result);
usage("reference, tdna, flank, reads file all need to exist") 
if any { ! defined $_ || ! -f $_ } ($reference_file, $tdna_file, $flank_file, $reads_file);

LOG("tdna_file: $tdna_file");
LOG("reference_file: $reference_file");
LOG("flank_file: $flank_file");
LOG("reads_file: $reads_file");

usage("need --output-directory (-d)") if ! $output_directory;

sub usage { pod2usage(scalar(@_) ? (-msg => shift()) : (), -verbose => 2, -noperldoc => 1); }


make_path $output_directory;

my $scaffold_file = do{
    my $tdna_file_basename  = basename($tdna_file, qw{.fa .fas .fasta});
    my $flank_file_basename = basename($flank_file, qw{.fa .fas .fasta});
    catdir $output_directory, "scaffold-$tdna_file_basename-$flank_file_basename.fasta";
};

my $output_file  = catfile $output_directory, basename($reads_file, qw{.fq .fastq}) . "_vs_" .  basename($scaffold_file, '.fasta') . ".sam";
#my $summary_file = catfile $output_directory, basename($output_file, ".sam") . ".summary.txt");

my $full_read_len = fastq_read_length($reads_file);
my $read_len      = @splice == 2 ? $splice[1] - $splice[0] + 1 : $full_read_len;

my $trim5 = @splice ? $splice[0] - 1              : 0;
my $trim3 = @splice ? $full_read_len - $splice[1] : 0;

$flank_prefix_length //= int($read_len * .75);
$tdna_prefix_length  //= int($read_len * .75);

LOG("read lengths: $full_read_len ($read_len used)");
LOG("read splice: @splice (which is trim5 $trim5 trim3 $trim3)");
LOG("flank_prefix_length: $flank_prefix_length");
LOG("tdna_prefix_length:  $tdna_prefix_length");

sub LOG { 
    my $msg = shift;
    say STDERR "LOG: $msg";
}

#######################################################################

my %scaffold;

# 0. Read reference

my $reference_fr = FastaReader->new(file => $reference_file, slurp => 1);
bowtie_build( file => $reference_file, version => 2);

# 1. Find where the flanks are from

LOG "aligning each flank to genome";

my $flank_fr = FastaReader->new(file => $flank_file);
my %flanks;

for my $s ($flank_fr->sequence_list()) {
    my $flank_seq = $flank_fr->get($s);
    if (length($flank_seq) < $flank_prefix_length){
        LOG "flank $s is too short, skipping";
    }
    else{
        my $alignments = Run::Bowtie2::where_is($reference_file, $flank_seq);
        for my $align (@$alignments) {
            LOG "$s: $align";
        }
        $flanks{$s} = [grep { $_->{mapped} } @$alignments];
    }
}

# 2. Get tDNA's left border's prefix, rc'd.

my $tdna_fr  = FastaReader->new(file => $tdna_file);
die "$tdna_file has more than one sequence?" if $tdna_fr->sequence_count != 1;

my $tdna_prefix = $tdna_fr->get($tdna_fr->first_sequence(), 1, $tdna_prefix_length, rc => 1);
LOG "tdna_prefix is $tdna_prefix";

# 3. Combine tdna prefix with flanks.

for my $f (keys %flanks) {
    my $flank_prefix = $flank_fr->get($f, 1, $flank_prefix_length);
    $scaffold{"tdna+$f"} = $tdna_prefix . $flank_prefix;

    # for each flank alignment, get the upstream region in the genome
    # and add to scaffold
    for my $l (@{$flanks{$f}}) {
        my $leftmost = $l->{leftmost};
        my $rightmost = $l->{rightmost};
        my $seqid = $l->{seqid};
        my $strand = $l->{reverse} ? "R" : "F";

        if ($l->{reverse}){
            my $upstream_start = $rightmost + 1;
            my $upstream_end = $rightmost + $tdna_prefix_length;

            $scaffold{"flank_upstream_${seqid}_${upstream_start}_${upstream_end}_${strand}+$f"} = 
            $reference_fr->get($seqid, $upstream_start, $upstream_end, rc => 1) . $flank_prefix;
        }
        else{
            my $upstream_start = $leftmost - $tdna_prefix_length;
            my $upstream_end = $leftmost - 1;

            $scaffold{"flank_upstream_${seqid}_${upstream_start}_${upstream_end}_${strand}+$f"} = 
            $reference_fr->get($seqid, $upstream_start, $upstream_end, rc => 1) . $flank_prefix;
        }
    }
}

# 4. Print scaffold and build index

LOG("producing scaffold file $scaffold_file");

open my $fh, '>', $scaffold_file;

for my $name (sort keys %scaffold) {
    LOG ">$name";
    LOG "$scaffold{$name}";
    say $fh ">$name";
    say $fh "$scaffold{$name}";
}

close $fh;

LOG("building bisulfite bowtie index for $scaffold_file");

my ($scaffold_bsrc_file) = bowtie_build( file => $scaffold_file, bs => 'c2t', rc => 1, version => 2,);

# 5. Convert reads 

my $cmd = "perl -S fastq2rcfasta.pl --c2t $reads_file | bowtie2 --norc -x $scaffold_bsrc_file -U - -f -S $output_file -5 $trim5 -3 $trim3 @ARGV";
LOG("running: $cmd");
system($cmd);

#######################################################################
# split the reads in smaller chunks (while converting to fasta/c2t)

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




=cut
