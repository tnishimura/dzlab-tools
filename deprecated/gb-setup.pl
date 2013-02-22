#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use YAML qw/LoadFile/;

my $config = shift // die "$0 conf.yaml";
my $staging_config = LoadFile($config)->{stagingconf};
my $gbrowse_config = LoadFile($config)->{gbrowseconf};
say $config;
say $staging_config;
say $gbrowse_config;

system("gb-prepare.pl $config");
system("gb-create-config.pl $staging_config > $gbrowse_config");
system("gb-load-mysql.pl -c $staging_config");

