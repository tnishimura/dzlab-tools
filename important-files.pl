#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use File::Find;
use File::Basename;
use File::Spec::Functions;
use Pod::Usage;
use Getopt::Long;

my $windows;
my $singlec;
my $inter;
my $divorce;
my $copy;
my $result = GetOptions (
    "copy|c=s" => \$copy,
    "single-c|s" => \$singlec,
    "windows|w"  => \$windows,
    "intermidiate|i" => \$inter,
    "divorce|d" => \$divorce,
);
if (!@ARGV || !$result){
    say "usage: ";
    say "$0 [-s] [-w] [-i] [-d] dir1 dir2 file1 file2 ... | rsync -avnrP --files-from=- . /somewhere";
    exit 1;
}

my %accum;
my @dirs = grep { -d } @ARGV;
my @files = grep { -f } @ARGV;

find( sub {
        # $File::Find::name - filename relative to pwd
        # $File::Find::dir  - dirname relative to pwd 
        # $_                - filename relative to $File::Find::dir

        my $dirbase = basename($File::Find::dir);

        # default extensions
        if (/\.log$/ || /\.table.txt$/ || /\.freq$/ || /\.ends$/ || 
            /\.avg$/ || /ratio\.txt$/ || /ratio\.txt\.no_coord_check$/ ||
            /log\.txt$/ ||
            # source code
            /\.pl$/ || /\.sh$/ 
        )
        { 
            $accum{$File::Find::name} = 1;
        }
        # if dir is windows is named 
        if ($windows && $dirbase =~ /^windows/){
            $accum{$File::Find::name} = 1;
        }
        # merged single-c files
        if ($singlec && $dirbase =~ /^single-c/ && /gff.merged$/){
            $accum{$File::Find::name} = 1;
        }
        # non-single-c/non-window big files
        if ($inter && $dirbase !~ /^windows/ && $dirbase !~ /^single-c/){ 
            if (/gff$/ || /c2t$/ || /g2a$/ || /fa$/ || /eland3$/ || /eland3\.post$/){
                $accum{$File::Find::name} = 1;
            }
        }

        if ($divorce && ( /\.7\.w50.*\.gff$/ || /\.3\.elfiltered$/ || /\.4\.gff$/)){
            $accum{$File::Find::name} = 1;
        }
    }, @dirs) if @dirs;


my $total_size = 0;
for my $file (@files, sort keys %accum) {
    my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size, $atime,$mtime,$ctime,$blksize,$blocks) = stat($file);
    $total_size += $size;
    say $file;
}

$total_size /= 1024*1024;
say STDERR "Total size: $total_size MB";

if (defined $copy && -d $copy){
    open my $rsync, '|-', "rsync -avWP --files-from=- . \"$copy\"";
    for my $file (@files, sort keys %accum) {
        say $rsync $file;
    }
    close $rsync;
}
