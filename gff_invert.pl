#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use List::Util qw/max min/;
use List::MoreUtils qw/all/;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use FastaReader;

END {close STDOUT}
$| = 1;

sub invert{
    my ($list, $start, $end) = @_;

    die "argument error" unless ref $list eq 'ARRAY';

    my @inverted;
    
    for my $element (@$list) {
        push @inverted, [$start, $element->[0]];
        $start = $element->[1];
    }
    push @inverted, [$start, $end];

    return \@inverted;
}

# input:
#   $list      - arrayref to list of elems to be merge
#   $get_start - sub to get start coord of $list elem
#   $get_end   - sub to get end coord of $list elem
# output:
#   [
#     [start1, end1, @list_of_elems_contributing_to_this_interval_1],
#     [start2, end2, @list_of_elems_contributing_to_this_interval_2],
#     ...
#   ]
sub bubble_merge{
    my ($list, $get_start, $get_end) = @_;

    die "argument error" unless all { ref $_ eq 'CODE' } ($get_start, $get_end);
    die "argument error" unless ref $list eq 'ARRAY';

    # interval version of cmp/<=>
    #   obj1 is [start, end, @elems],
    #   obj2 is [start, end, $elem]
    # return:
    #   -1 if $obj1 is behind $obj2
    #   0  if $obj1 is overlapping $obj2
    #   1  if $obj1 is after $obj2
    my $cmp_by_interval = sub {
        my ($obj1, $obj2) = @_;
        my ($start1, $end1, $start2, $end2) = ($obj1->[0], $obj1->[1], $obj2->[0], $obj2->[1]);

        if ($end1 < $start2){
            return -1;
        }
        elsif ($end2 < $start1){
            return 1;
        }
        else{
            return 0;
        }
    };

    # merge interval, concat element list
    #   in_progess is [start, end, @elems]
    #   current is    [start, end, $elems]
    #   output is     [start, end, @elems]
    my $interval_merge = sub{
        my ($prev, $current) = @_;
        my ($start1, $end1, @accum) = @$prev;
        my ($start2, $end2, $newelem) = @$current;

        return [min($start1, $start2), max($end1, $end2), @accum, $newelem];
    };


    # convert lists to [start, end, element], sort by start
    my @sorted = sort { 
        $a->[0] <=> $b->[0]
    } 
    map { 
        [$get_start->($_), $get_end->($_), $_] 
    } @$list;

    my @merged;

    # in_progess is [start, end, @elems]
    my $in_progress = shift @sorted; # first one.

    for my $current (@sorted) {
        my $cmp = $cmp_by_interval->($in_progress, $current);
        if ($cmp == -1){ # if $current is complete past $in_progress
            push @merged, $in_progress;
            $in_progress = $current;
        }
        elsif ($cmp == 1){ # $current behind $in_progress shouldn't happen
            die "wtf? not sorted?";
        }
        else{ # overlapping-- merge to $in_progress
            $in_progress = $interval_merge->($in_progress, $current);
        }
    }
    push @merged, $in_progress; # last one
    return \@merged;
}

sub gff_get_start{ $_[0]->start(); }
sub gff_get_end{ $_[0]->end(); }

sub run{
    use Pod::Usage;
    use Getopt::Long;
    
    my $result = GetOptions (
        "reference|r=s" => \my $reference,
        "feature|f=s"   => \(my $feature = "inversion"),
    );

    if (!$result || ! $reference){
        say "usage: gff_invert.pl ...";
        exit 1;
    }

    my $fasta = FastaReader->new(file => $reference);
    my $parser = GFF::Parser->new(file => \*ARGV, normalize => 0);

    my %slurp;

    while (defined(my $gff = $parser->next())){
        my $seq = $gff->sequence();
        next if ! $seq;
        push @{$slurp{$seq}}, $gff;
    }

    my $id = "0000000";
    for my $seq (sort keys %slurp) {
        #my ($list, $cmp_by_start, $cmp_by_interval, $interval_merge) = @_;
        my $merged = bubble_merge($slurp{$seq}, \&gff_get_start, \&gff_get_end);
        my $inverted = invert($merged, 0, $fasta->get_length($seq));
        for my $block (@$inverted) {
            my ($start, $end) = ($block->[0]+1, $block->[1]-1);
            if ($start <= $end){
                say join "\t", $seq, q{.}, $feature, $start, $end, qw{. + .}, "ID=I" . ++$id;
            }
        }
    }
}

run() unless caller();
