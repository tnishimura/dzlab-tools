#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More qw(no_plan);
use Test::Exception;

use DZUtil qw/c2t rc_c2t g2a rc_g2a/;
use Run::BowtieBuild;
use TestUtils;
use FastaReader;

my $ref = setup_reference();
my ($filename, $index_prefix, @ebwt) = bowtie_build(file => $ref, force => 1);
my ($c2t_filename, $c2t_index_prefix, @c2t_ebwt) = bowtie_build(file => $ref, bs => 'c2t', force => 1);
my ($g2a_filename, $g2a_index_prefix, @g2a_ebwt) = bowtie_build(file => $ref, bs => 'g2a', force => 1);

# make sure it's RC'ing and bs'ing correctly
my $fr = FastaReader->new(file => $filename);
my $c2t_fr = FastaReader->new(file => $c2t_filename);
my $g2a_fr = FastaReader->new(file => $g2a_filename);

my @expected_seqs = sort map { $_, "RC_$_" } $fr->sequence_list();
is_deeply([sort $c2t_fr->sequence_list()], \@expected_seqs, "c2t correct sequences");
is_deeply([sort $g2a_fr->sequence_list()], \@expected_seqs, "g2a correct sequences");

for my $seqid ($fr->sequence_list()) {
    my $original = $fr->get($seqid);
    is(c2t($original),    $c2t_fr->get($seqid));
    is(rc_c2t($original), $c2t_fr->get("RC_" . $seqid));
    is(g2a($original),    $g2a_fr->get($seqid));
    is(rc_g2a($original), $g2a_fr->get("RC_" . $seqid));
}

for my $ebwt (@ebwt, @c2t_ebwt, @g2a_ebwt) {
    ok(-f $ebwt && -s $ebwt);
}

