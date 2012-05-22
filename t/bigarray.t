#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More qw(no_plan);
use Test::Exception;
use BigArray;
use PDL;
use List::MoreUtils qw//;

#######################################################################
# randomized assignment 

SKIP: {
    #skip "not now";

    my @success;
    my $N = 10_000;
    my $ntest = 20;

    for (1 .. $ntest) {
        my $ba = BigArray->new(size => $N, base => 0);
        my $aref = [ (0) x $N ];
        for (1 .. 20 * $N) {
            my $index = int(rand($N));
            my $value = int(rand(100));

            $ba->push_pair($index, $value);
            $ba->push_increment($index);

            $aref->[$index] += 1;
            $aref->[$index] += $value;
        }
        #$ba->commit_increment();
        #$ba->commit_pair();
        $ba->commit();
        push @success, sum_square($ba->get_pdl(), pdl($aref)) == 0;
    }
    all_ok(\@success, "randomized assignments");
}

sub sum_square{
    my ($x, $y) = @_;
    return (($x - $y) ** 2)->sum();
}

sub all_ok{
    my ($aref, $msg) = @_;
    ok(List::MoreUtils::all(sub{ $_ }, @$aref), $msg);
}

#######################################################################
# commit_increment

{
    my $insert_count = 100_000;
    my $ba_size      = 10_000;
    my $ba_bufsize   = int($ba_size / 10);
    my $ba           = BigArray->new(size => $ba_size, base => 0, buffer_size => $ba_bufsize);

    my @residence_ok;
    for my $num_pushed (0 .. $insert_count - 1) {
        my $residence = $ba->get_buffer_residence;
        if ($num_pushed % $ba_bufsize == 0){
            push @residence_ok, $residence == 0;
            #is($residence, 0, "residence $num_pushed");
        }
        $ba->push_increment(int(rand($ba_size)));
    }
    all_ok(\@residence_ok, "committing appropriately");
}

#######################################################################
# create_collection

my %specs = (
    chr1 => 100,
    chr2 => 200,
    chr3 => 123,
);

my $collection = BigArray::create_collection({}, \%specs);
while (my ($k,$v) = each %specs) {
    is($collection->{$k}{size}, $specs{$k}, "create_collection size test");
}
