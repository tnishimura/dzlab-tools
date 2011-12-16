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

sub bubble_merge{
    my ($list, $get_start, $get_end) = @_;

    die "argument error" unless all { ref $_ eq 'CODE' } ($get_start, $get_end);
    die "argument error" unless ref $list eq 'ARRAY';

    my @sorted = 
    sort { 
        $a->[0] <=> $b->[0]
    } 
    map {
        [ $get_start->($_), $get_end->($_) ]
    } 
    @$list;

    my @merged;
    my $in_progress = shift @sorted;

    for my $current (@sorted) {
        my $cmp = interval_cmp($in_progress, $current);
        if ($cmp == -1){ 
            push @merged, $in_progress;
            $in_progress = $current;
        }
        elsif ($cmp == 1){
            die "wtf? not sorted?";
        }
        else{
            $in_progress = interval_merge($in_progress, $current);
        }
    }
    push @merged, $in_progress;
    return \@merged;
}

sub interval_merge{
    return [min($_[0]->[0], $_[1]->[0]), max($_[0]->[1], $_[1]->[1])];
}

sub interval_cmp{
    my ($object1, $object2) = @_;
    my ($start1, $end1, $start2, $end2) = (@$object1, @$object2);

    if ($end1 < $start2){
        return -1;
    }
    elsif ($end2 < $start1){
        return 1;
    }
    else{
        return 0;
    }
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
        if (! exists $slurp{$seq}){
            $slurp{$seq} = [];
        }
        push @{$slurp{$seq}}, $gff;
    }

    my $id = "0000000";
    for my $seq (sort keys %slurp) {
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
