#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More;
use Test::Exception;
use YAML qw/Load Dump LoadFile DumpFile/;
use Sam::Parser;
use Readonly;

# basics
{
    my $parser = Sam::Parser->new(
        file => "t/data/sam-test.sam", 
        skip_unmapped => 0,
    );

    my $count = 0;
    my $mapped = 0;
    while (defined(my $sam = $parser->next())){
        $count++;
        $mapped++ if $sam->mapped;
    }
    is($count, 100, "correct read count");
    is($mapped, 84, "correct mapped read count");
}

# RC_
{
    my $parser = Sam::Parser->new(
        file => "t/data/sam-test.sam", 
        skip_unmapped => 0,
        convert_rc => 1
    );

    my $rc_seq_count = 0;
    while (defined(my $sam = $parser->next())){
        $rc_seq_count ++ if $sam->seqid =~ /^RC_/;
    }
    is($rc_seq_count, 0, "no RC_ found");
}

done_testing();
