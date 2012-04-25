#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More qw(no_plan);
use Test::Exception;
use TestUtils;
use FastaUtil;
use FastaReader;
use DZUtil qw/reverse_complement c2t g2a/;

my $ref = setup_reference;
say $ref;

my $rc_file    = "$ref.rc";
my $c2t_file   = "$ref.c2t_test";
my $g2a_file   = "$ref.g2a_test";
my $c2trc_file = "$ref.rc.c2t_test";
my $g2arc_file = "$ref.rc.g2a_test";

rc_fasta_on_disk($ref, $rc_file);
bs_fasta_on_disk('c2t',$ref, $c2t_file);
bs_fasta_on_disk('g2a',$ref, $g2a_file);
bsrc_fasta_on_disk('c2t',$ref, $c2trc_file);
bsrc_fasta_on_disk('g2a',$ref, $g2arc_file);

my $fr       = FastaReader->new(file => $ref);
my $rc_fr    = FastaReader->new(file => $rc_file);
my $c2t_fr   = FastaReader->new(file => $c2t_file);
my $g2a_fr   = FastaReader->new(file => $g2a_file);
my $c2trc_fr = FastaReader->new(file => $c2trc_file);
my $g2arc_fr = FastaReader->new(file => $g2arc_file);

is(scalar($rc_fr->sequence_count()), 14);
is(scalar($c2t_fr->sequence_count()), 7);
is(scalar($g2a_fr->sequence_count()), 7);
is(scalar($c2trc_fr->sequence_count()), 14);
is(scalar($g2arc_fr->sequence_count()), 14);

my @expected_seqs = sort qw/chr1 chr2 chr3 chr4 chr5 chrm chrc/;
my @expected_rc_seqs = map { "RC_$_" } @expected_seqs;
my @expected_all_seqs = sort @expected_seqs, @expected_rc_seqs;

is_deeply(\@expected_all_seqs, [sort $rc_fr->sequence_list]);
is_deeply(\@expected_seqs,     [sort $c2t_fr->sequence_list]);
is_deeply(\@expected_seqs,     [sort $g2a_fr->sequence_list]);
is_deeply(\@expected_all_seqs, [sort $c2trc_fr->sequence_list]);
is_deeply(\@expected_all_seqs, [sort $g2arc_fr->sequence_list]);

for my $seq (@expected_seqs) {
    is($fr->get($seq),      $rc_fr->get($seq), "rc $seq");
    is(c2t($fr->get($seq)), $c2t_fr->get($seq), "c2t $seq");
    is(g2a($fr->get($seq)), $g2a_fr->get($seq), "g2a $seq");
    is(c2t($fr->get($seq)), $c2trc_fr->get($seq), "c2trc $seq");
    is(g2a($fr->get($seq)), $g2arc_fr->get($seq), "g2arc $seq");

    my $rcseq = "RC_$seq";

    is($fr->get($seq, undef, undef, rc => 1),      $rc_fr->get($rcseq), "rc $rcseq");
    is(c2t($fr->get($seq, undef, undef, rc => 1)), $c2trc_fr->get($rcseq), "c2trc $rcseq");
    is(g2a($fr->get($seq, undef, undef, rc => 1)), $g2arc_fr->get($rcseq), "g2arc $rcseq");
}
