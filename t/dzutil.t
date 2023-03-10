#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/../lib";
use DZUtil qw/memofile fastq_read_length overlap common_suffix common_prefix chext split_names/;
use Test::Exception;


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


ok(common_prefix( "11234_", "1123_", "112_",) eq "112", 'common_prefix');
ok(common_prefix( "", "1123_", "112_",) eq "", 'common_prefix');
ok(common_prefix( 123, 12223, "112_",) eq 1, 'common_prefix');

ok(common_suffix( "11234_", "1123_", "112_",) eq "_", 'common_suffix');
ok(common_suffix( "", "1123_", "112_",) eq "", 'common_suffix');
ok(common_suffix( 123, 12223, "11asdf223",) eq '23', 'common_suffix');

ok(overlap([1,40], [39, 45]) == 2, 'overlap');

is(100, fastq_read_length("t/data/fastq_read_length/fastq_read_length-test.fastq"), "fastq_read_length plain");
is(100, fastq_read_length("t/data/fastq_read_length/fastq_read_length-test.fastq.gz"), "fastq_read_length gzip");
is(100, fastq_read_length("t/data/fastq_read_length/fastq_read_length-test.fastq.bz2"), "fastq_read_length bzip2");
is(100, fastq_read_length("t/data/fastq_read_length/dir-gzip"), "fastq_read_length gzip dir");
is(100, fastq_read_length("t/data/fastq_read_length/dir-bz2"), "fastq_read_length bzip2 dir");
dies_ok { fastq_read_length "t/data/fastq_read_length/dir-gzip-uneven" } "fastq_read_length gzip dir uneven dies";

### memofile 
is(memofile("/home/neko/meow.txt", '/tmp'), "/tmp/home,neko,meow.txt", "memofile");
