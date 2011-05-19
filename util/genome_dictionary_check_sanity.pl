#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use YAML qw/LoadFile DumpFile/;
use Log::Log4perl qw/:easy/;
Log::Log4perl->easy_init( { 
    level    => $DEBUG,
    #file     => ">run.log",
    layout   => '%d{HH:mm:ss} %p> (%L) %M - %m%n',
} );

my $counter=0;
my $badcounter=0;
my $logger = get_logger();
for my $file (@ARGV) {
    my $config = LoadFile($file);

    my %align = %{$config->{alignment}};

    while (my ($ch,$alignments) = each %align) {
        for my $quad (@$alignments) {
            $counter++;
            my ($x1,$x2, $y1,$y2) = @$quad;
            if ($x2-$x1 == $y2-$y1){
                #$logger->info("$ch @$quad");
            } else{
                $logger->info("$ch @$quad = " . ($x2-$x1 - ($y2-$y1)));
                $badcounter++;
                #$logger->info($x2-$x1 - ($y2-$y1));
            }
        }
    }
}


$logger->info("$badcounter / $counter alignments didn't match in size");
