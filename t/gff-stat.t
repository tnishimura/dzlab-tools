#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More qw(no_plan);
use GFF::Statistics;
use List::Util qw/sum/;

#######################################################################
# Utility tests

sub mean{
    return @_ ? sum(@_)/scalar(@_) : 'na';
}
sub median{
    if (@_ == 0){
        return 'na';
    }
    my @s = sort { $a <=> $b } @_;
    return $s[int(@s/2)];
}
my $tol = .0001;

for (1 .. 10){
    my @test = map { int(rand(100)) } (1..1000);
    my %testh = GFF::Statistics::tohist(@test);
    ok(mean(@test) - GFF::Statistics::histmean(\%testh) < $tol, 'histmean');
    is(median(@test), GFF::Statistics::histmedian(\%testh), 'histmedian');
    is(sum(@test), GFF::Statistics::sumhists(\%testh), 'sumhists');
    is(2*sum(@test), GFF::Statistics::sumhists(\%testh, \%testh), 'sumhists with 2');
}

for (1 .. 10) {
    my @test = map { rand(100) } (1..1000);
    my $avg = GFF::Statistics::make_averager();
    for (@test){
        $avg->($_);
    }
    ok($avg->() - mean(@test) < $tol, "make_averager");
}

#######################################################################
# getstats

sub safemeth{
    my ($c, $t) = @_;
    return ($c + $t > 0) ? ($c / ($c + $t)) : 'na';
}

sub naive_getstats{
    my $file = shift;
    my $p = GFF::Parser->new(file => $file);
    my (@chr_ct, @mit_ct, @nuc_ct, @chr_methyl, @mit_methyl, @nuc_methyl);
    my %ct = (
        chr => {c => 0, t => 0},
        mit => {c => 0, t => 0},
        nuc => {c => 0, t => 0},
    );
    GFF:
    while (defined(my $gff = $p->next())){
        given ($gff->sequence){
            my ($c, $t) = ($gff->get_column('c'), $gff->get_column('t'));
            next GFF if ! defined $c || ! defined $t || $t + $c <= 0;
            my $ct = $c + $t;
            my $methyl = $c / $ct;
            when (/chrc/i){
                push @chr_ct, $ct;
                push @chr_methyl, $methyl;
                $ct{chr}{c}+=$c;
                $ct{chr}{t}+=$t;
            }
            when (/chrm/i){
                push @mit_ct, $ct;
                push @mit_methyl, $methyl;
                $ct{mit}{c}+=$c;
                $ct{mit}{t}+=$t;
            }
            when (/chr\d+/i){
                push @nuc_ct, $ct;
                push @nuc_methyl, $methyl;
                $ct{nuc}{c}+=$c;
                $ct{nuc}{t}+=$t;
            }
        }
    }
    return {
        chr_ct_median   => median(@chr_ct),
        mit_ct_median   => median(@mit_ct),
        nuc_ct_median   => median(@nuc_ct),
        chr_ct_mean     => mean(@chr_ct),
        mit_ct_mean     => mean(@mit_ct),
        nuc_ct_mean     => mean(@nuc_ct),
        chr_methyl_mean => mean(@chr_methyl),
        mit_methyl_mean => mean(@mit_methyl),
        nuc_methyl_mean => mean(@nuc_methyl),
        chr_methyl_total => safemeth($ct{chr}{c}, $ct{chr}{t}),
        mit_methyl_total => safemeth($ct{mit}{c}, $ct{mit}{t}),
        nuc_methyl_total => safemeth($ct{nuc}{c}, $ct{nuc}{t}),
        coverage        => sum(@chr_ct, @mit_ct, @nuc_ct),
    };
}


my $test_cg = "t/gff-stat-cg-test-data.gff";
my $test_chh = "t/gff-stat-chh-test-data.gff";

sub compare{
    my ($x, $y) = @_;
    if ($x eq 'na' && $y eq 'na'){ return 1; }
    elsif ($x eq 'na' && $y ne 'na'){ return 0; }
    elsif ($x ne 'na' && $y eq 'na'){ return 0; }
    else{
        return $x - $y < $tol;
    }
}
for my $file ($test_cg, $test_chh) {
    my $stats = getstats($file);
    my $naive = naive_getstats($file);
    for my $rowname (@GFF::Statistics::rownames) {
        ok(compare($stats->{$rowname}, $naive->{$rowname}), "naive vs smart: $rowname");
    }
}

#say Dumper getstats($test_cg);
#say Dumper naive_getstats($test_cg);
#say Dumper getstats($test_chh);
