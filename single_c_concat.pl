#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/lib";
use autodie;
use Config::General qw(ParseConfig);
use Data::Dumper;
use feature 'say';
use File::Basename;
use File::Find;
use File::Spec::Functions;
use DZUtil qw/localize common_prefix common_suffix timestamp/;

die "usage: $0 contextual-single-c-files.gff ..." unless @ARGV;

my $prefix = common_prefix(@ARGV);
my $suffix = common_suffix(@ARGV);
$prefix =~ s/\.$|\_$//;
$suffix =~ s/^\.|\^_//;
my $concat  = $prefix . ".all." . $suffix;

open my $outfh, '>', $concat;
for my $f (@ARGV) {
    open my $infh, '<', $f;
    while (defined(my $line = <$infh>)){
        chomp $line;
        say $outfh $line;
    }
    close $infh;
}
close $outfh;
