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
my $result = GetOptions (
    "single-c|s" => \$singlec,
    "windows|w"  => \$windows,
);
if (!@ARGV || !$result){
    say "usage: ";
    say "$0 [-s] [-w] dir1 dir2 file1 file2 ... | rsync -avnrP --files-from=- . /somewhere";
    exit 1;
}

my %accum;
my @dirs = grep { -d } @ARGV;
my @files = grep { -f } @ARGV;

find( sub {
        # $File::Find::name - filename relative to pwd
        # $File::Find::dir  - dirname relative to pwd 
        # $_                - filename relative to $File::Find::dir

        # default extensions
        if (/log$/ || /table$/ || /freq$/ || /ends$/ || /avg$/ || /ratio\.txt$/){ 
            $accum{$File::Find::name} = 1;
        }
        # if dir is windows is named 
        if (($windows) && basename($File::Find::dir) =~ /^windows$/){
            $accum{$File::Find::name} = 1;
        }
        if (($singlec) && basename($File::Find::dir) =~ /^single-c$/ && /gff.merged$/){
            $accum{$File::Find::name} = 1;
        }
    }, @dirs) if @dirs;

for my $file (@files, sort keys %accum) {
    say $file;
}

