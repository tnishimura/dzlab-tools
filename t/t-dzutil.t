#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/../lib";
use DZUtil qw/chext split_names/;


use Test::More qw(no_plan);

is(chext("/etc/passwd.txt", "bat"), "/etc/passwd.bat", "chext 1 abs");
is(chext("home/asdf/passwd.gff", "eland"), "home/asdf/passwd.eland", "chext 2 rel");
is(chext("home/asdf/passwd", "eland"), "home/asdf/passwd.eland", "chext 3 without initial extension");
is_deeply(
    [split_names("/etc/passwd.txt",qw/a b c d ee/)],
    [qw{
    /etc/passwd-a.txt
    /etc/passwd-b.txt
    /etc/passwd-c.txt
    /etc/passwd-d.txt
    /etc/passwd-ee.txt
    }],
    "split_names 1",
);

is_deeply(
    [split_names("/etc/passwd",qw/a b c d ee/)],
    [qw{
    /etc/passwd-a
    /etc/passwd-b
    /etc/passwd-c
    /etc/passwd-d
    /etc/passwd-ee
    }],
    "split_names 2",
);

