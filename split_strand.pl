#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Cwd qw/getcwd/;

END {close STDOUT}
$| = 1;

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    #"longname|l=isf" => \my $var,
    "col|c=i" => \my $col,
);
if (!$result || ! defined $col){
    say "usage: bowtie_split_strand.pl ...";
    exit 1;
}

my $bowtie = shift || die "usage: $0 file";

my $basename = $bowtie;
my $ext = '';
if ($basename =~ /(\.\w+)$/){
    $ext = $1;
}

open my $fh, '<', $bowtie;
open my $plus, '>', "$basename.plus$ext";
open my $minus, '>', "$basename.minus$ext";

while (defined(my $line = <$fh>)){
    chomp $line;
    my @F = split /\t/, $line;
    if ($F[$col] eq '+'){
        say $plus $line;
    }
    elsif ($F[$col] eq '-'){
        say $minus $line;
    }
}
close $fh;
close $minus;
close $plus;

