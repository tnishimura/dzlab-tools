#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use IO::File;
use IO::Handle;

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "help"    => \my $help,
);

my ($left_file, $right_file) = @ARGV;
my $paired = defined $right_file;

my %handles;
$handles{left} = IO::File->new($left_file);
$handles{right} = IO::File->new($right_file) if $paired;

while (my ($side,$fh) = each %handles) {
    my $strand = $paired && $side eq 'left' ? '+' : 
                 $paired && $side eq 'right' ? '-' : 
                 '.';

    my $seq;
    my $feature;
    while (defined(my $line = <$fh>)){
        chomp $line;
        if ($line =~ /^track/){
            ($feature) = $line =~ /name="(\w+)"/;
        }
        elsif ($line =~ /^variableStep/){
            ($seq) = $line =~ /chrom=(\w+)/;
            say $seq;
        }
        elsif (! defined $seq){
            die "data before track info?"
        }
        else{
            my ($pos, $score) = split /\s+/, $line;
            say join("\t", $seq, '.', $feature // '.', $pos, $pos, $score, $strand, '.', '.');
        }
    }
}
__DATA__
track	type=wiggle_0	name="Col0_rep1 minus strand"	visibility=full	color=20,150,20	altColor=150,20,20	windowingFunction=mean
variableStep	chrom=chr1
110	1
116	1
117	0.285714
125	0.714286
126	0.571429
129	0.428571
139	0.428571
141	0.428571
