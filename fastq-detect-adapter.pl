#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;
use List::Util qw/first max min shuffle sum/;

my $result = GetOptions (
    "adapter-length|l=i" => \(my $adapter_length),
    "top|t=i" => \(my $top),
    "max-reads|m=i" => \(my $max_reads),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || ! $adapter_length);  

my %adapter_count;
my %adapter_count_per_read;
my $total = 0;
my $read_count = 0;

while (<ARGV>){
    my $read = <ARGV>;
    <ARGV>;
    <ARGV>;
    chomp ($read);

    my %found_in_this_read;

    for my $i (0 .. (length($read) - $adapter_length)) {
        my $x = substr($read, $i, $adapter_length);
        $total += 1;
        $adapter_count{$x}++;
        if (! $found_in_this_read{$x}){
            $adapter_count_per_read{$x}++;
            $found_in_this_read{$x} = 1;
        }
    }
    $. % 250_000 == 0 and warn $.;

    if ($max_reads < $read_count++){
        last;
    }
}

# my @sorted = sort { $b->[1] <=> $a->[1] } map { [$_, $adapter_count{$_} / $read_count] } keys %adapter_count;
my @sorted          = sort { $b->[1] <=> $a->[1] } map { [$_, $adapter_count{$_} / $total] } keys %adapter_count;
my @sorted_per_read = sort { $b->[1] <=> $a->[1] } map { [$_, $adapter_count_per_read{$_} / $read_count] } keys %adapter_count_per_read;

say "rank\tsequence\tshare of sequences of length $adapter_length";
for my $level (0 .. min($top - 1, $#sorted)) {
    say "$level\t$sorted[$level][0]\t$sorted[$level][1]";
}
say "";

say "rank\tsequence\tnumber of reads with sequence";
for my $level (0 .. min($top - 1, $#sorted_per_read)) {
    say "$level\t$sorted_per_read[$level][0]\t$sorted_per_read[$level][1]";
}
