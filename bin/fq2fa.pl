#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/../lib";
use DZUtil qw/open_maybe_compressed/;

END {close STDOUT}
$| = 1;

for my $file (@ARGV) {
    my $fh = open_maybe_compressed($file);

    while (defined(my $line = <$fh>)){
        chomp $line;
        if ($line =~ s/^@//){
            my ($first, $second) = split ' ', $line;

            # most common
            if ($first =~ m{#0/[12]$}){
                print ">$first\n";
            }
            # new solexa format, 6/20/2011
            elsif (defined $second && $second =~ /^([12])/){
                print ">$first#/$1\n";
            }
            # handle older reads with no read id's 
            else{
                print ">$first#/1\n";
            }
        }
        $_ = <$fh>; print;
        <$fh>; <$fh>;
    }
    close $fh;
}

