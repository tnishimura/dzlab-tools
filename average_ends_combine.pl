#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use DZUtil qw/common_suffix common_prefix/;
use Tie::IxHash;
use File::Basename qw/basename/;

use Getopt::Long;

my $result = GetOptions (
    "strip-suffix|s" => \(my $strip_suffix),
    "strip-prefix|p" => \(my $strip_prefix),
);
usage() if (!$result);  

END {close STDOUT}
$| = 1;

usage() if ! @ARGV;
my $prefix = common_prefix map { basename($_) } @ARGV;
my $suffix = common_suffix map { basename($_) } @ARGV;

my %accum;
tie %accum, "Tie::IxHash";

for my $file (@ARGV) {
    my $key = basename($file);
    $key =~ s/^\Q$prefix\E// if $strip_prefix;
    $key =~ s/\Q$suffix\E$// if $strip_suffix;

    if (exists $accum{key}){
        die "duplicate column names";
    }
    $accum{$key} = slurp_avg($file)
}

my @lens = map { scalar @{$accum{$_}} } keys %accum;
my $len = shift @lens;
if (grep { $len ne $_ } @lens){
    warn "all files not same length";
}

say join "\t", keys %accum;

for my $i (0 .. $len - 1) {
    say join "\t", map { $_->[$i] // '.' } values %accum;
}

#######################################################################
# util

sub slurp_avg{
    my $file = shift ;
    my @accum;
    open my $fh, '<', $file;
    <$fh>;
    while (defined(my $line = <$fh>)){
        chomp $line;
        push @accum, [split /\t/, $line]->[1];
    }

    close $fh;
    return \@accum;
}

sub usage{
    say "$0 [-p] [-s] *.ends.avg > combined.ends.avg";
    exit 1;
}
