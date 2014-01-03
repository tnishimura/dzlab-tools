#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use File::Basename qw/basename/;
use Getopt::Long qw/:config no_ignore_case/;
use List::MoreUtils qw/any/;
use Pod::Usage;
use IO::All;
use File::Spec::Functions qw/catfile/;
use Parallel::ForkManager;

use FindBin;
use lib "$FindBin::Bin/../lib";

use FastqReader::Convert;
use Run::BowtieBuild;
use Run::Bowtie2;

my $result = GetOptions (
    "reference-file|r=s"   => \(my $reference_file),
    "flank-file|f=s"       => \(my $flank_file),
    "reads-file|R=s"       => \(my $reads_file),
    "insert-seq|s=s"       => \(my $insert_seq),
    "insert-pos|p=i"       => \(my $insert_pos),
    "armor-size|a=i"       => \(my $armor_size = 50),
    "output-directory|o=s" => \(my $output_directory),
    "parallel|t=i" => \(my $parallel),
);

my $pm = Parallel::ForkManager->new($parallel);

usage() unless $reference_file && $flank_file && $reads_file && $insert_seq && $insert_pos && $output_directory;

#######################################################################
# setup output directory and files

my $outdir = io($output_directory)->mkpath();

my $flank_basename     = basename($flank_file, qw{.fa .fasta});
my $reference_basename = basename($reference_file, qw{.fa .fasta});
my @upstream_range   = ($insert_pos - $armor_size + 1, $insert_pos);
my @downstream_range = ($insert_pos + 1, $insert_pos + $armor_size - 1);
my @across_range     = ($insert_pos - $armor_size, $insert_pos + $armor_size);

my $upstream_with_flank   = catfile($outdir,"${reference_basename}_$upstream_range[0]_$upstream_range[1].$flank_basename.fasta");
my $downstream_with_flank = catfile($outdir,"$flank_basename.${reference_basename}_$downstream_range[0]_$downstream_range[1].fasta");
my $across                = catfile($outdir,"${reference_basename}_$across_range[0]_$across_range[1].fasta");

my $reads_bsrc_file = $reads_file . ".c2t";

my $upstream_output   = catfile $outdir, basename($reads_file) . '-vs-' . basename($upstream_with_flank, '.fasta') . '.sam';
my $downstream_output = catfile $outdir, basename($reads_file) . '-vs-' . basename($downstream_with_flank, '.fasta') . '.sam';
my $across_output     = catfile $outdir, basename($reads_file) . '-vs-' . basename($across, '.fasta') . '.sam';

#######################################################################

my $flank_fr     = FastaReader->new(file => $flank_file);
my $reference_fr = FastaReader->new(file => $reference_file);

# say Dumper $insert_seq, $reference_fr, $reference_file, $reference_fr->get($insert_seq, @upstream_range);

my $flank_seq = $flank_fr->get($flank_fr->first_sequence());

io($upstream_with_flank  )->print(">up\n" .   $reference_fr->get($insert_seq, @upstream_range)  . $flank_seq );
io($downstream_with_flank)->print(">down\n" . $flank_seq . $reference_fr->get($insert_seq, @downstream_range));
io($across)->print(">across\n" . $reference_fr->get($insert_seq, @across_range));

my ($upstream_with_flank_c2t)   = bowtie_build( file => $upstream_with_flank,   bs => 'c2t', rc => 1, version => 2);
my ($downstream_with_flank_c2t) = bowtie_build( file => $downstream_with_flank, bs => 'c2t', rc => 1, version => 2);
my ($across_c2t) = bowtie_build( file => $across, bs => 'c2t', rc => 1, version => 2);

#######################################################################

# LOG("additional bowtie2 options: @ARGV");

# convert reads c2t
if (! -f $reads_bsrc_file || -s $reads_bsrc_file < .25 * -s $reads_file){
    fastq_convert(to_fasta => 1, methyl => 'c2t', in => $reads_file, out => $reads_bsrc_file);
}

if ($pm->start == 0){
    bowtie2_raw(
        [ 
            -U => $reads_bsrc_file,
            -x => $upstream_with_flank_c2t,
            -S => $upstream_output,
            '-f',
        ],
        [
            convert_rc => 1,
            skip_unmapped => 1,
        ],
    );
    $pm->finish; 
}

if ($pm->start == 0){
    bowtie2_raw(
        [ 
            -U => $reads_bsrc_file,
            -x => $downstream_with_flank_c2t,
            -S => $downstream_output,
            '-f',
        ],
        [
            convert_rc => 1,
            skip_unmapped => 1,
        ],
    );
    $pm->finish; 
}

if ($pm->start == 0){
    bowtie2_raw(
        [ 
            -U => $reads_bsrc_file,
            -x => $across_c2t,
            -S => $across_output,
            '-f',
        ],
        [
            convert_rc => 1,
            skip_unmapped => 1,
        ],
    );
    $pm->finish; 
}

$pm->wait_all_children;

#######################################################################
# aux function

sub LOG { 
    my $msg = shift;
    io("=")->println("LOG: $msg\n");
    # io($log_file)->appendln("LOG: $msg\n");
}
sub usage { 
    pod2usage(scalar(@_) ? (-msg => shift()) : (), -verbose => 2, -noperldoc => 1); 
}


