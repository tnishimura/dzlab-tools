#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use Fastq;

die "$0 file.fastq id" unless @ARGV == 2;

my $f = Fastq->new(file => $ARGV[0]);

say $f->get($ARGV[1]);



