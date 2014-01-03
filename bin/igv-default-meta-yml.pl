#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use IO::All;
use Cwd qw/getcwd/;
use File::Basename qw/basename dirname/;
use File::Path qw/make_path remove_tree/;
use File::Spec::Functions qw/rel2abs canonpath catdir catfile updir/;
use File::Copy;
use File::Temp qw/tempdir/;


my @files = @ARGV;

use YAML qw/Load Dump LoadFile DumpFile/;

die "meta.yml already exists? delete manually if you really want to overwrite" if -f 'meta.yml';
DumpFile('meta.yml', {
        name =>  'test',
        annotations =>[
            map {
                {
                    name => basename($_) =~ s/\.\w+$//r,
                    path => $_,
                }
            } @files
        ]
    });

