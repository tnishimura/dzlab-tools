#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More qw(no_plan);
use Test::Exception;
use CHI::Memo;

use File::Temp qw/tempdir/;
my $tempdir = tempdir(CLEANUP => 1);

$CHI::Memo::root_dir = $tempdir;
CHI::Memo::clear_cache();

my $calc_counter = 0;

for (0 .. 10) {
    chimemo '/etc/passwd', 'size', sub { $calc_counter++; -s $_ };
}

is($calc_counter, 1, "calc_counter");

my $hash_payload = {a => 1, b => 2, c => [1,2,3] };
my $array_payload = [1,3,4, { x => 1, y => 654} ];
my $scalar_payload = 42;
my @list_payload = ($hash_payload, $array_payload, $scalar_payload);

chimemo '/etc/passwd', 'hash_payload', sub { $calc_counter++; $hash_payload };
chimemo '/etc/passwd', 'hash_payload', sub { $calc_counter++; $hash_payload };
is($calc_counter, 2, "calc_counter");

chimemo '/etc/passwd', 'array_payload', sub { $calc_counter++; $array_payload };
chimemo '/etc/passwd', 'array_payload', sub { $calc_counter++; $array_payload };
is($calc_counter, 3, "calc_counter");

chimemo '/etc/passwd', 'list_payload', sub { $calc_counter++; @list_payload };
chimemo '/etc/passwd', 'list_payload', sub { $calc_counter++; @list_payload };
is($calc_counter, 4, "calc_counter");

chimemo '/etc/passwd', 'scalar_payload', sub { $calc_counter++; $scalar_payload };
chimemo '/etc/passwd', 'scalar_payload', sub { $calc_counter++; $scalar_payload };
is($calc_counter, 5, "calc_counter");

is_deeply(
    chimemo('/etc/passwd', 'hash_payload', sub { $calc_counter++; $hash_payload }),
    $hash_payload, "hash payload",
);

is_deeply(
    chimemo('/etc/passwd', 'array_payload', sub { $calc_counter++; $array_payload }),
    $array_payload, "array payload",
);

is_deeply(
    chimemo('/etc/passwd', 'scalar_payload', sub { $calc_counter++; $scalar_payload }),
    $scalar_payload, "scalar payload",
);

is_deeply(
    [chimemo('/etc/passwd', 'list_payload', sub { $calc_counter++; @list_payload })],
    \@list_payload, "list payload",
);

is($calc_counter, 5, "calc_counter");
