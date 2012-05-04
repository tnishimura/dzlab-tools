#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
END {close STDOUT}
$| = 1;

use Pod::Usage;
use Getopt::Long;
use File::Path qw/make_path/;
use File::Spec::Functions qw/rel2abs canonpath catdir catfile updir/;
use Tie::File;
use File::Basename qw/basename dirname/;

use FindBin;
use lib "$FindBin::Bin/lib";
use BowtieParser;
use FastaReader;
use Digest::MD5::Util;

my $result = GetOptions (
    "tmpdir|d=s"    => \my $tmpdir,
    "reference|r=s" => \my $ref,
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || ! $tmpdir || ! $ref);

#######################################################################

make_path($tmpdir);

my $fr = FastaReader->new(file => $ref, slurp => 0);
my $bp = BowtieParser->new(file => \*ARGV);

my $checksum_file = catfile($tmpdir, "MD5SUM");

#######################################################################
# read and split

my @filenames;   # filename

if (-f $checksum_file && defined(my $confirmed_files = md5_confirm($checksum_file))){
    @filenames = @$confirmed_files;
}
else{
    my %filehandles; # seq => fh

    while ( 
        defined(my $left_bowtie = $bp->next()) and 
        defined(my $right_bowtie = $bp->next())
    ){
        my %left = process($left_bowtie);
        my %right = process($right_bowtie);

        my $seq = $left{seq};

        if (! exists $filehandles{$seq}){
            my $filename = catfile($tmpdir, $seq);
            open my $fh, '>', $filename;
            $filehandles{$seq} = $fh;
            push @filenames, $filename;
        }

        say {$filehandles{$seq}} join "\t", $left{start}, $right{end};
    }
    close $_ for values %filehandles;

    for my $file (@filenames) {
        system("sort -k1,1 -k2,2n -i $file -o $file");
    }

    md5_create_checksum($checksum_file, @filenames);
}

#######################################################################

for my $file (@filenames) {
    my $seq = filename2seq($file);

    my $ptr = 0;
    my $length = $fr->get_length($seq);

    tie my @lines, 'Tie::File', $file;
    my $numlines = @lines;

    for my $chr_position (1 .. $length) {
        #say "=== $seq $chr_position ===";
        my @overlappers;
        my @overlapper_midpoints;
        my $lookahead = 0;
        FINDOVERLAPS:
        while ($ptr + $lookahead  < $numlines){
            my ($start, $end) = split /\t/, $lines[$ptr + $lookahead];
            if ($end < $chr_position){ # position is already past this read, increment
                $ptr++;
            }
            elsif ($chr_position < $start){
                last FINDOVERLAPS;
            }
            else{
                push @overlappers, [$start, $end];
                push @overlapper_midpoints, ($end + $start) / 2;
                $lookahead++;
            }
        }

        if (@overlapper_midpoints){
            #say Dumper \@overlappers;
            say join "\t", 
            $seq, qw/. ./, $chr_position, $chr_position, 
            normal_mad(\@overlapper_midpoints), qw/. ./, "n=" . scalar(@overlapper_midpoints);
        }
    }
}

#######################################################################
# utilities

sub seq2filename{
    my ($tmpdir, $seq) = @_;
    return catfile($tmpdir, $seq);
}

sub filename2seq{
    my ($file) = @_;
    basename $file;
}

sub normal_mad{
    my $aref = shift;
    my $med = median($aref);
    return 1.4826 * median([map { abs($_ - $med) } @$aref]);
}

sub median{
    my $aref = shift;
    if (@$aref == 1){
        return $aref->[0];
    }
    elsif (@$aref == 0){
        die "empty aref passed to median";
    }
    my @sorted = sort { $a <=> $b } @$aref;
    my $numelems = scalar @sorted;
    if ($numelems % 2 == 0){
        return ($sorted[($numelems - 1) / 2] + $sorted[($numelems - 1) / 2 + 1])/2;
    }
    else{
        return $sorted[($numelems - 1) / 2];
    }
}

sub process{
    my $bowtie = shift;
    return 
        seq   => $bowtie->[2],
        start => $bowtie->[3],
        end   => $bowtie->[3] + length($bowtie->[4]) - 1;
}
