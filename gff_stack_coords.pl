#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use File::Basename qw/basename/;
use Scalar::Util qw/looks_like_number/;
use Pod::Usage;
END {close STDOUT}
$| = 1;

use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "gap|g=i" => \(my $gap = 5000),
);

if (! $result || (! @ARGV && -t STDIN)){
    pod2usage(-verbose => 2, -noperldoc => 1);
}
my $col1 = @ARGV == 1 ? basename($ARGV[0], '.gff') : 'stacked';

my $p = GFF::Parser->new(file => \*ARGV);
my $last;

my $offset = 0;
while (defined(my $gff = $p->next())){
    my ($seq, $orig_start, $orig_end, $score) = $gff->slice(1,4,5,6);
    $score //= '.';

    # invalid line
    if (! defined($seq) || ! looks_like_number($orig_start) || ! looks_like_number($orig_end)){
        next;
    }

    if (! $last){
        $last = $gff;
    }
    elsif ( $last->sequence() ne $seq ){
        $offset = $last->end() - $orig_start + $gap + 1;
    }

    my $start = $orig_start + $offset;
    my $end   = $orig_end + $offset;

    say join "\t", $col1, qw/. ./, $start, $end, $score, qw/. ./, "seq=$seq;start=$orig_start;end=$orig_end";
    $last = $gff;
}

=head1 NAME

gff_stack_coords.pl - Convert chromosome coordinates. written as one-shot for
yvonne so she could feed multi-chromosome GFF file to autocorrelation. Assumes
inputs are sorted.  Place --gap between chromosomes (default 5000);

=head1 SYNOPSIS

Usage examples:

 gff_stack_coords.pl file.gff > output.gff
 gff_stack_coords.pl < file.gff > output.gff
 some_other_command.pl | gff_stack_coords.pl > output.gff

With custom gap:

 gff_stack_coords.pl --gap 10000 file.gff > output.gff
 gff_stack_coords.pl -g 1000 < file.gff > output.gff

=cut

