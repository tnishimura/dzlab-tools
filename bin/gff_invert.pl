#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use List::Util qw/max min/;
use List::MoreUtils qw/all/;
use FindBin;
use lib "$FindBin::Bin/../lib";
use GFF::Parser;
use FastaReader;
use IntervalMerge;

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
        my $merged = interval_merge($slurp{$seq}, \&gff_get_start, \&gff_get_end);
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
