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

use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use GFF::Statistics qw/gff_detect_width/;
use Conjure;

my $result = GetOptions (
    "reference-genomee|r=s" => \(my $reference_genomee),
    "igv-tools-jar|j=s" => \(my $igv_tools_jar = "$FindBin::Bin/share/igvtools.jar"),
    # "score-from-ct|ct" => \(my $score_from_ct),
    "tmp-dir|d=s" => \(my $tmp_dir),
    "no-feature-in-filename|nff" => \(my $no_feature_in_filename),
    "parallel|p=i" => \(my $parallel=0),
);
pod2usage(-verbose => 2, -noperldoc => 1) 
if (!$result || ! $igv_tools_jar || ! $reference_genomee || ! -f $igv_tools_jar || ! -f $reference_genomee);  

my $pm = Parallel::ForkManager->new($parallel);

# GFF high-density single-c/windows files
for my $methylgff_file (@ARGV){
    $pm->start and next;
    my $width = gff_detect_width($methylgff_file);
    my $infh = IO::File->new($methylgff_file);

    # write each feature/seqid combo to seperate tmpfile, combine later
    my %tmpfiles; # {methyl|coverage}{feature}{seqid}
    while (defined(my $gffline = $infh->getline)){
        chomp $gffline;
        my ($seq, undef, $feature, $start, $end, $score, undef, undef, $attr) = split /\t/, $gffline;

        use Scalar::Util qw/looks_like_number/;
        # say $score if ! looks_like_number($score);
        if (looks_like_number($score) and all { $_ ne '.' } $seq, $feature, $start, $end, $attr){
            if (! exists $tmpfiles{methyl}{$feature}{$seq}){
                $tmpfiles{methyl}{$feature}{$seq} = File::Temp->new(DIR => $tmp_dir, UNLINK => 1);
                $tmpfiles{coverage}{$feature}{$seq} = File::Temp->new(DIR => $tmp_dir, UNLINK => 1);
                my $mtmp = $tmpfiles{methyl}{$feature}{$seq}->filename;
                my $ctmp = $tmpfiles{coverage}{$feature}{$seq}->filename;
                warn "creating tmp files for $feature $seq ($mtmp, $ctmp)";
                # variableStep chrom=chr1  span=50
                # 0	0.083500
                # 50	0.139000
                # 100	0.119200
                # 150	0.142900
                # 200	0.035300
                # 250	0.026800
                # 300	0.041700
                # $tmpfiles{methyl}{$feature}{$seq}->print("variableStep chrom=$seq  span=$width\n");
                # $tmpfiles{coverage}{$feature}{$seq}->print("variableStep chrom=$seq  span=$width\n");
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
    }
    # combine each type/feature combination
    for my $type (qw/methyl coverage/) {
        my %feature2seq2fh = %{$tmpfiles{$type}};
        while (my ($feature,$seq2fh_ref) = each %feature2seq2fh) {
            my $outfile = $no_feature_in_filename ?  $methylgff_file =~ s/\.gff$/.$type.wig/r : $methylgff_file =~ s/\.gff$/.$feature-$type.wig/r;
            my $outfh = IO::File->new($outfile, 'w');

            for my $seq (sort keys %$seq2fh_ref) {
                my $tmpfh = $tmpfiles{$type}{$feature}{$seq};
                my $tmpfile = $tmpfh->filename;
                $tmpfh->close;

                # sort
                warn "sorting $tmpfile";
                conjure program => "sort -k1,1n -i $tmpfile -o $tmpfile";

                # append
                my $sortedfh = IO::File->new($tmpfile);
                $outfh->print("variableStep chrom=$seq span=$width\n");
                while (defined(my $line = <$sortedfh>)){
                    $outfh->print($line);
                }
            }
            close $outfh;
            my $tdf = $outfile =~ s/\.wig$/\.tdf/r;
            conjure(program => "java -jar $igv_tools_jar toTDF $outfile $tdf $reference_genomee",
                on_stderr => sub { say },
                on_stdout => sub { say },
            );
        }
    }
    $pm->finish; 
}
$pm->wait_all_children;

=head1 gff2tdf.pl 

Usage examples:

 gff2tdf.pl [-p #threads] -r reference.fas input.gff

=cut

