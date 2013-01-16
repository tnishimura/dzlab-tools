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

my $original_rc_count = 0;

# basics
{
    my $parser = Sam::Parser->new(
        file          => "t/data/sam-test.sam",
        skip_unmapped => 0,
        convert_rc    => 0,
    );

    my $count = 0;
    my $mapped = 0;
    my $reversed = 0;
    while (defined(my $sam = $parser->next())){
        $count++;
        $mapped++ if $sam->mapped;
        $reversed++ if $sam->is_reverse;
        $original_rc_count++ if $sam->seqid =~ /^RC_/;
    }
    is($count, 100, "correct read count");
    is($mapped, 84, "correct mapped read count");
    is($reversed, 0, "should be no reverse strand");
}

# RC_
{
    my $parser = Sam::Parser->new(
        file => "t/data/sam-test.sam", 
        skip_unmapped => 0,
        convert_rc => 1
    );

    my $rc_count = 0;
    my $reverse_strand_count = 0;
    my $reverse_indicator = 0;
    while (defined(my $sam = $parser->next())){
        $rc_count++             if $sam->seqid =~ /^RC_/;
        $reverse_strand_count++ if $sam->is_reverse;
        $reverse_indicator++    if $sam->mapped && $sam->readid =~ /:\-:/;
    }
    is($rc_count, 0, "no RC_ found");
    is($original_rc_count, $reverse_strand_count, "RC_'s were properly turned to - strand");
    ok($reverse_indicator >= $reverse_strand_count, "RC_ count is appropriate as indicated by readid");
    # >= here instead of == b/c one of the reads aligned elsewhere
}

done_testing();
