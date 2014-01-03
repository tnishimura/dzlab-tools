#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Pod::Usage;
use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/../lib";
use GFF::Parser;
use IntervalMerge qw/interval_merge interval_deinterlace/;
use GFF::Slurp qw/gff_slurp_index/;

END {close STDOUT}
$| = 1;

my $result = GetOptions (
    "ignore-feature|i=s" => \my @ignore,
    "output|o=s"         => \(my $output = '-'),
);
pod2usage(-verbose => 1) if (!$result || !@ARGV);

if ($output ne '-'){
    open my $fh, '>', $output;
    select $fh;
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

{
    my $id = "0000000";
    sub id{
        return "ID=M" . ++$id;
    }
}

sub gff_get_start { $_[0]->start() }
sub gff_get_end { $_[0]->end() }

for my $seq (sort keys %seq2gff) {
    my $gffs = $seq2gff{$seq};
    say STDERR "$seq interval_merge start. number of elements: ", scalar(@$gffs);
    my $merged = interval_merge($gffs, \&gff_get_start, \&gff_get_end);

    say STDERR "$seq interval_merge done. dumping";
    for my $island (@$merged) {
        my ($start, $end, @elems) = @$island;

        my @features = all_features(@elems);
        my $new_feature = @features == 1 ? $features[0] : 'mixed';

        for my $piece (@{interval_deinterlace(\@elems, \&gff_get_start, \&gff_get_end)}) {
            say join "\t",
            $seq,
            '.',
            $new_feature,
            $piece->[0], 
            $piece->[1],
            '.',
            '+',
            '.',
            id();
        }
    }

    say STDERR "dumping done";
}

=head1 NAME

gff_merge_annotation.pl - Merge overlapping windows in an annotation GFF file.

=head1 SYNOPSIS

Basic usage:

 gff_merge_annotation.pl -o deinterlaced.gff yourannotation.gff 

Some annotations include 'chromosome' features which will screw everything up
(b/c it overlaps with everything), so you'll want to ignore it:

 gff_merge_annotation.pl -o merged.gff -i chromosome TAIR8_gmod.gff 

=head1 OPTIONS

=over

=item  -i <feature> | --ignore-feature <feature>

If given, lines with this feature will be ignore.  For example, '-i gene' will
ignore all genes in the file.  You can specifiy multiple times if desired.


=back

=cut
