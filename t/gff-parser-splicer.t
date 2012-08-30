#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Carp;

use FindBin;
use lib "$FindBin::Bin/../lib";

use GFF::Parser;
use GFF::Parser::Splicer;
use TestUtils;
use Test::GFF;
use Test::More;

my @default_colnames = qw/sequence start end score c t n/;

my $testgff = test_gff(1_000);
my $splicer = GFF::Parser::Splicer->new(file => $testgff, columns => \@default_colnames);
my $parser  = GFF::Parser->new(file => $testgff, slip => 1);

my %splicer_stats = my %parser_stats = (c => 0, t => 0, n => 0, score => 0);

while(1){
    my $splicer_gff = $splicer->next();
    my $parser_gff = $parser->next();

    if (!defined $splicer_gff && !defined $parser_gff){
        last;
    }
    elsif (defined $splicer_gff xor defined $parser_gff){
        croak "not giving same number of lines?";
    }

    my ($sequence, $start, $end, $score, $c, $t, $n) = @$splicer_gff;
    $splicer_stats{c} += $c;
    $splicer_stats{t} += $t;
    $splicer_stats{n} += $n;
    $splicer_stats{score} += $score;

    $parser_stats{c} += $parser_gff->get_column('c');
    $parser_stats{t} += $parser_gff->get_column('t');
    $parser_stats{n} += $parser_gff->get_column('n');
    $parser_stats{score} += $parser_gff->get_column('score');
}

is_deeply(\%splicer_stats, \%parser_stats, "splicer and parser give same results");

done_testing();
