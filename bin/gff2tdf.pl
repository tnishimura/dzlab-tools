#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use IO::File;
use File::Temp;
use List::MoreUtils qw/all/;
use Parallel::ForkManager;
use Pod::Usage;
use Getopt::Long;
use Scalar::Util qw/looks_like_number/;

use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use FastaReader;
use GFF::Statistics qw/gff_detect_width/;

my $result = GetOptions (
    "reference-genomee|r=s"      => \(my $reference_genomee),
    "igv-tools-jar|j=s"          => \(my $igv_tools_jar = "$FindBin::Bin/share/igvtools.jar"),
    "tmp-dir|d=s"                => \(my $tmp_dir),
    "no-feature-in-filename|nff" => \(my $no_feature_in_filename),
    "parallel|p=i"               => \(my $parallel=0),
    "score-only|s"               => \(my $score_only),
);
pod2usage(-verbose => 2, -noperldoc => 1) 
if (!$result || ! $igv_tools_jar || ! $reference_genomee || ! -f $igv_tools_jar || ! -f $reference_genomee);  

my $ref = FastaReader->new(file => $reference_genomee);
my $pm = Parallel::ForkManager->new($parallel);

# GFF high-density single-c/windows files
for my $gff_file (@ARGV){
    $pm->start and next;
    my $width = gff_detect_width($gff_file);
    my $infh = IO::File->new($gff_file);

    # write each feature/seqid combo to seperate tmpfile, combine later
    my %tmpfiles; # {score|methyl|coverage}{feature}{seqid}
    while (defined(my $gffline = $infh->getline)){
        chomp $gffline;
        my @fields = split /\t/, $gffline;
        next if @fields != 9;

        my ($seq, undef, $feature, $start, $end, $score, undef, undef, $attr) =  @fields;
        $seq = $ref->get_original_name($seq);

        # methylation is the default
        if (! $score_only and all { $_ ne '.' } $seq, $feature, $start, $end, $attr){
            if (! exists $tmpfiles{methyl}{$feature}{$seq}){
                $tmpfiles{methyl}{$feature}{$seq} = File::Temp->new(DIR => $tmp_dir, UNLINK => 1);
                $tmpfiles{coverage}{$feature}{$seq} = File::Temp->new(DIR => $tmp_dir, UNLINK => 1);
                my $mtmp = $tmpfiles{methyl}{$feature}{$seq}->filename;
                my $ctmp = $tmpfiles{coverage}{$feature}{$seq}->filename;
                warn "creating tmp files for $feature $seq ($mtmp, $ctmp)";
            }
            my ($c, $t);
            if ($attr=~/c=(\d+)/){ $c = $1; }
            if ($attr=~/t=(\d+)/){ $t = $1; }
            if (defined $c && defined $t){
                my $coverage = $c + $t;
                next if $coverage == 0;
                my $methylation = $c / $coverage;
                $tmpfiles{methyl}{$feature}{$seq}->print("$start\t$methylation\n");
                $tmpfiles{coverage}{$feature}{$seq}->print("$start\t$coverage\n");
            }
        }

        # --score-only is an addition to support non-methyl plain gff files with the score in col 6
        # feature does not need to be non-'.' for score_only
        if ($score_only and looks_like_number($score) and all { $_ ne '.' } $seq, $start, $end){
            if (! exists $tmpfiles{score}{$feature}{$seq}){
                $tmpfiles{score}{$feature}{$seq} = File::Temp->new(DIR => $tmp_dir, UNLINK => 1);
                warn "creating tmp score file for $feature $seq";
            }
            $tmpfiles{score}{$feature}{$seq}->print("$start\t$score\n");
        }
    }
    # combine each type/feature combination
    for my $type (qw/score methyl coverage/) {
        next if ! exists $tmpfiles{$type};

        my %feature2seq2fh = %{$tmpfiles{$type}};
        while (my ($feature,$seq2fh_ref) = each %feature2seq2fh) {
            my $wigfile = $gff_file;
            if ($no_feature_in_filename || $feature eq '.'){
                $wigfile =~ s/\.gff$/.$type.wig/; 
            }
            else{
                $wigfile =~ s/\.gff$/.$feature-$type.wig/;
            }
            my $outfh = IO::File->new($wigfile, 'w');

            # sort the tmp files, and append them to $wigfile
            # yes, all sequences can be in the same file
            for my $seq (sort keys %$seq2fh_ref) {
                my $tmpfh = $tmpfiles{$type}{$feature}{$seq};
                my $tmpfile = $tmpfh->filename;
                $tmpfh->close;

                # sort
                warn "sorting $tmpfile";
                system "sort -k1,1n -i $tmpfile -o $tmpfile";

                # append
                my $sortedfh = IO::File->new($tmpfile);
                $outfh->print("variableStep chrom=$seq span=$width\n");
                while (defined(my $line = <$sortedfh>)){
                    $outfh->print($line);
                }
                close $sortedfh;
            }
            close $outfh;

            my $tdf = $wigfile;
            $tdf =~ s/\.wig$/\.tdf/;
            system "java -jar $igv_tools_jar toTDF $wigfile $tdf $reference_genomee";
        }
    }
    $pm->finish; 
}
$pm->wait_all_children;

=head1 gff2tdf.pl 

This script converts a gff file into a TDF file suitable for loading into IGV.
It requires that java is installed and available in your PATH.  It can handle
methylation files (which is the default), which case it produces two tdf files
(for methylation and coverage), or plain GFF files, in which case a single tdf
file is created from the scores column (column 6).

This creates input.CG-methyl.tdf, input.CG-coverage.tdf, input.CG-methyl.tdf, etc..

 gff2tdf.pl -r reference.fas input.gff

This creates input.score.tdf:

 gff2tdf.pl -s -r reference.fas input.gff

In both cases, you can add a "-p #" argument to run in parallel.

You should not redirect the output (don't use ">").

=cut

# wig files look like this:

# variableStep chrom=chr1  span=50
# 0	0.083500
# 50	0.139000
# 100	0.119200
# 150	0.142900
# 200	0.035300
# 250	0.026800
# 300	0.041700
# variableStep chrom=chr2  span=50
# 0	0.083500
# 50	0.139000
# 100	0.119200
# 150	0.142900
# 200	0.035300
# 250	0.026800
# 300	0.041700
