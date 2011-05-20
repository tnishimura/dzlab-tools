#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Test::More qw(no_plan);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Launch;
use Log::Log4perl qw/get_logger/;
use File::Path qw/rmtree/;

my $conf=q/
    log4perl.logger          = DEBUG, Print
    log4perl.appender.Print                          = Log::Log4perl::Appender::Screen
    log4perl.appender.Print.layout                   = PatternLayout
    log4perl.appender.Print.layout.ConversionPattern = m %d{HH:mm:ss} %p> (%L) %M - %m%n
/;
Log::Log4perl::init( \$conf );

ok(launch('sleep 1') );

