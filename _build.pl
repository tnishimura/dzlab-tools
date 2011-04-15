#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
if (@ARGV != 1){
    print "usage: $0 tag-name";
    exit 255;
}

#rm -rf help
# pod2projdocs -o help -l . -except doc/index.pod -except '^tmp'
system("perl -wlpe 's/__VERSION__/$ARGV[0]/' installer.nsi.in > installer.nsi");
system("perl -wlpe 's/__VERSION__/$ARGV[0]/' dzlab-check.pl.in > dzlab-check.pl");
system("makensis installer.nsi");
