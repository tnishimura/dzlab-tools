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

my $result = GetOptions (
    "keep-going|k" => \my $keepgoing,
);
if (!$result or (!@ARGV and -t STDIN)){
    say "Detect windows line ending and non-ascii characters";
    say "usage: non-ascii.pl [-k | --keep-going] file";
    exit 1;
}

my %notified;
my %warned;
while (defined(my $line = <ARGV>)){
    if (! exists $notified{$ARGV}){
        say STDERR "Processing $ARGV";
        $notified{$ARGV}++;
    }

    my $ending = $line =~ tr/\r\n//d;
    if ($ending == 2){
        $warned{$ARGV} = 1;
    }
    if (0 < $line =~ tr/\t\x20-\x7e//c){
        if ($keepgoing){
            if (! exists $warned{$ARGV}){
                warn "$ARGV\n$.\n$line";
                $warned{$ARGV}++;
            }
        }
        else{
            die "$ARGV $.\n$line";
        }
    }

    my $numcols = scalar split /\t/, $line;
    my $is_blank = $line =~ /^\s*$/;
    my $is_comment = $line =~ /^\s*#$/;

    if ($numcols != 9 and ! $is_blank and ! $is_comment){
        die "$ARGV $.\n$line";
    }
}

if (keys %warned){
    say STDERR join "\n", "following files had windows line endings:", 
    sort keys %warned;
}
