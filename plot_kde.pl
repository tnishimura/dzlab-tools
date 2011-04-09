#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Scalar::Util qw(looks_like_number );

use Pod::Usage;
use Getopt::Long;

my $opt_min;
my $opt_max;
my $opt_buckets = 50;
my $output = '-';
my $fh;

my $result = GetOptions (
    "buckets|b=i" => \$opt_buckets,
    "max=f" => \$opt_max,
    "min=f" => \$opt_min,
    "output|o=s" => \$output,
);
pod2usage(-verbose => 1) if (!$result || ! @ARGV);  

if ($output eq '-'){
    $fh = \*STDOUT;
} else {
    open $fh, '>', $output;
}

my @accum;

while (defined(my $line = <ARGV>)){
    chomp $line;
    $line =~ tr/\r\n//d;
    if (looks_like_number $line){
        push @accum, $line;
    } else {
        warn "Non-numeric line \"$line\" on line $. of $ARGV";
    }
}
if (! @accum) { die "file empty?"; }

@accum = sort { $a <=> $b } @accum;
my $num = scalar @accum;
my $min = $opt_min // $accum[0];
my $max = $opt_max // $accum[-1];
my $bucket_size = ($max - $min ) / $opt_buckets;
my @buckets = ();

for my $bucki (0 .. $opt_buckets - 1) {
    my $left = $min + $bucket_size * $bucki;
    my $right = $min + $bucket_size * ($bucki + 1);
    my $total = 0;
    while (@accum && $accum[0] < $right){
        $total++;
        shift @accum;
    }
    printf $fh "%f\t%f\t%f\n", $left,$right, ($total / $num / $bucket_size);
}

if ($output ne '-'){
    close $fh;
}


=head1 NAME

plot_kde.pl - Feed this script a list of numbers to produce a histogram.

=head1 SYNOPSIS

Usage examples:

 plot_kde.pl -b 100 -o histogram.txt data.txt

=head1 OPTIONS

=over

=item  --max <num>

Manually specify max.  Default to the max actually found in list.

=item  --min <num>

Manually specify min.  Default to the min actually found in list.

=item  -b <num> | --buckets <num>

Number of buckets in histogram. Default 50.

=item  -o <file> | --output <file>

Output file, default to screen (STDOUT).

=back

=cut

