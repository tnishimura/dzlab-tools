#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;

use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use IntervalMerge qw/interval_merge/;
use GFF::Slurp qw/gff_slurp_index/;

END {close STDOUT}
$| = 1;

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "ignore-feature|i=s" => \my @ignore,
);
if (!$result){
    say "usage: gff_merge_annotation.pl ...";
    exit 1;
}

sub all_features{
    map { $_->feature() } @_;
}

my %seq2gff;
my $parser = GFF::Parser->new(file => \*ARGV, normalize => 0, launder => 1);

say STDERR "slurping.";
while (defined(my $gff = $parser->next())){
    if (! defined $gff->feature() || ! grep { $_ eq $gff->feature } @ignore){
        push @{$seq2gff{$gff->sequence}}, $gff;
    }
}

say STDERR "done slurping.";

my $id = "0000000";
for my $seq (sort keys %seq2gff) {
    my $gffs = $seq2gff{$seq};
    say STDERR "$seq interval_merge start. number of elements: ", scalar(@$gffs);
    my $merged = interval_merge(
        $gffs, 
        sub { $_[0]->start() },
        sub { $_[0]->end() },
    );

    say STDERR "$seq interval_merge done. dumping";
    for my $island (@$merged) {
        my ($start, $end, @elems) = @$island;
        my @features = all_features(@elems);
        my $new_feature = @features == 1 ? $features[0] : 'mixed';
        say join "\t",
        $seq,
        '.',
        $new_feature,
        $start, 
        $end,
        '.',
        '+',
        '.',
        "ID=M" . ++$id;
    }
    say STDERR "dumping done";
}
