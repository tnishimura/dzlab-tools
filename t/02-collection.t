#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Tree::Range::Collection;
use Test::More qw(no_plan);

my $src = Tree::Range::Collection->new();

$src->add('seq1',10,20,"a");
$src->add('seq1',15,30,"b");
$src->add('seq1',30,50,"c");
$src->add('seq1',60,40,"d");
$src->add('seq1',80,70,"e");
$src->add('seq1',64,81,"f");

$src->add('seq2',10,20,"a");
$src->add('seq2',17,30,"b");
$src->add('seq2',30,50,"c");
$src->add('seq2',60,40,"d");
$src->add('seq2',80,70,"e");
$src->add('seq2',64,82,"f");
$src->finalize();
say $src->info;

is_deeply(
    [sort {$a->{item} cmp $b->{item}} $src->search_overlap('seq1',35,45)], 
    [sort {$a->{item} cmp $b->{item}} {item => "c", overlap => 11}, {item => "d", overlap => 6}], 
    "overlap 1"
);
is_deeply(
    [sort {$a->{item} cmp $b->{item}} $src->search_overlap('seq2',82,82)], 
    [sort {$a->{item} cmp $b->{item}} {item => "f", overlap => 1}], 
    "overlap 1"
);

ok(!$src->search('seq3',1,2), 'nonexistance');
