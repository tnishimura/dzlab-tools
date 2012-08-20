#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use File::Basename qw/basename/;
use Scalar::Util qw/looks_like_number/;
use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "gap|g=i" => \(my $gap = 5000),
);

if (! $result || (! @ARGV && -t STDIN)){
    pod2usage(-verbose => 2, -noperldoc => 1);
}
my $col1 = @ARGV == 1 ? basename($ARGV[0], '.gff') : 'stacked';

my $header_line = <ARGV>;
$header_line =~ tr/\n\r//d;
my @header = split /\t/, $header_line;
my $column_count = scalar @header;

say join "\t", @header, "abs_coord";

my $last;

my $offset = 0;
while (defined(my $line = <ARGV>)){
    $line =~ tr/\n\r//d;
    my ($seq, $orig_coord, @rest) = split /\t/, $line, $column_count;

    if (! defined($seq) || ! looks_like_number($orig_coord) ){
        die "line $.: invalid seq/coord? is this really a gff mashup file";
    }

    if (! $last){
        $last = [$seq, $orig_coord, @rest];
    }
    elsif ( $last->[0] ne $seq ){
        $offset = $last->[1] - $orig_coord + $gap + 1;
    }

    my $coord = $orig_coord + $offset;

    say join "\t", $seq, $orig_coord, @rest, $coord;
    $last = [$seq, $orig_coord, @rest];
}

=head1 NAME

gff_mashup_add_abs_coord.pl - Adds "chromosome coordinates". written as
one-shot for yvonne so she could feed multi-chromosome mashup file to
autocorrelation. Input must have seq and coord as first two columns (as it
should be for gff_mashup.pl output).  Place --gap between chromosomes (default
5000);

=head1 SYNOPSIS

Usage examples:

 gff_mashup_add_abs_coord.pl file.gff > output.gff
 gff_mashup_add_abs_coord.pl < file.gff > output.gff
 some_other_command.pl | gff_mashup_add_abs_coord.pl > output.gff

With custom gap:

 gff_mashup_add_abs_coord.pl --gap 10000 file.gff > output.gff
 gff_mashup_add_abs_coord.pl -g 1000 < file.gff > output.gff

=cut

