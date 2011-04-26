#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use File::Copy;
use FindBin;
use Pod::Usage;
use Getopt::Long;

my $bindir = $FindBin::Bin;
my $outdir = $bindir;
my $tag;

my $result = GetOptions (
    "outdir=s"  => \$outdir,
    "tag=s"     => \$tag,
);

if (!$result || !$tag){
    say "usage: $0 --tag tagname --outdir some/where/outdir";
    exit 1;
}

#rm -rf help
# pod2projdocs -o help -l . -except doc/index.pod -except '^tmp'
system("perl -wlpe 's/__VERSION__/$tag/' installer.nsi.in > installer.nsi");
system("perl -wlpe 's/__VERSION__/$tag/' dzlab-check.pl.in > dzlab-check.pl");
system("makensis installer.nsi");
unlink "installer.nsi";
unlink "dzlab-check.pl";

if ($outdir ne $bindir){
    move "dzlab-tools-$tag.exe", $outdir;
}

