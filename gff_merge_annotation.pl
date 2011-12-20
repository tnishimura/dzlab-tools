#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Pod::Usage;
use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use IntervalMerge qw/interval_merge/;
use GFF::Slurp qw/gff_slurp_index/;

END {close STDOUT}
$| = 1;

my $result = GetOptions (
    "ignore-feature|i=s" => \my @ignore,
    "locus-tag|t=s" => \my @locus_tags,
    "constituents|c" => \my $constituents,
);
pod2usage(-verbose => 1) if (!$result);

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
        my @elems = @_;
        if (! @locus_tags){
            return "ID=M" . ++$id;
        }
        else{
            my %tag_values;
            my $other;
            ELEM:
            for my $e (@elems) {
                for my $tag (@locus_tags) {
                    my $val = $e->get_column($tag);
                    if (defined $val){
                        for my $split_val (split /,/, $val) {
                            $tag_values{$split_val} = 1;
                        }
                        next ELEM;
                    }
                }
                $other = 1;
            }
            #die Dumper \@elems, \%tag_values;
            return "ID=" . join ",", (sort keys %tag_values), $other ? ('other') : ();
        }
    }

    sub constituents{
        my @elems = @_;
        my %features = map { $_ => 1 } all_features(@elems);
        return "Constituents=" . join ",", keys %features;
    }
}

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
        id(@elems) . ($constituents ? ";" . constituents(@elems) : ''),
    }

    say STDERR "dumping done";
}

=head1 NAME

gff_merge_annotation.pl - Merge overlapping windows in an annotation GFF file.

=head1 SYNOPSIS

Basic usage:

 gff_merge_annotation.pl yourannotation.gff > TAIR8_gmod.merged.gff

Some annotations include 'chromosome' features which will screw everything up
(b/c it overlaps with everything), so you'll want to ignore it:

 gff_merge_annotation.pl -i chromosome TAIR8_gmod.gff > TAIR8_gmod.merged.gff

If you want to gather the ID and Parent fields from entries, do (Note multiple
-t options):

 gff_merge_annotation.pl -i chromosome -t ID -t Parent TAIR8_gmod.gff > TAIR8_gmod.merged.gff

=head1 OPTIONS

=over

=item  -t <tag> | --locus-tag <tag>

If given, output GFF lines will have an 'ID=' field with collected ID's from
each constituent lines.  You can specifiy multiple times if desired.

=item  -i <varname> | --ignore-feature <varname>

If given, lines with this feature will be ignore.  For example, '-i gene' will
ignore all genes in the file.  You can specifiy multiple times if desired.

=item  -c | --constituents 

If given, output GFF lines will have an 'Constituent=' field with collected
features.

=back

=cut
