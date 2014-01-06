#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Pod::Usage;
use Getopt::Long;
END {close STDOUT}
$| = 1;

use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use FastaReader;

my $result = GetOptions (
    "plus"   => \my $plus,
    "minus"  => \my $minus,
    "dot"    => \my $dot,
    "prefix" => \(my $prefix = "RC_"),
);
my ($gff, $reference) = @ARGV;
if (!$result 
    || ! defined $gff 
    || ! defined $reference 
    || ! -f $gff 
    || ! -f $reference 
    || (1 < grep { defined $_ } ($plus, $minus, $dot))){
    pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]);
}

{
    my $default = $plus ? '+' : $minus ? '-' : $dot ? '.' : undef;
    sub strand{
        my $strand = shift;
        if ($default)             { return $default; }
        elsif (! defined $strand) { return '.'; }
        elsif ($strand eq '+')    { return '-'; }
        elsif ($strand eq '-')    { return '+'; }
        else                      { return '.'; }
    }
}

my $fr = FastaReader->new(file => $reference);
my $p = GFF::Parser->new(file => $gff, normalize => 0);

sub convert{
    my ($start, $end, $len) = @_;
    return ($len - $end + 1, $len - $start + 1);
}

while (defined(my $gff = $p->next())){
    my $seq = $gff->sequence();

    if (defined $seq && $seq =~ s/^$prefix//i){
        my $len = $fr->get_length($seq);
        say join "\t",
        $seq, 
        $gff->source()  // '.',
        $gff->feature() // '.',
        convert($gff->start(), $gff->end(), $len),
        $gff->score()            // '.',
        strand($gff->strand())   // '.',
        $gff->frame()            // '.',
        $gff->attribute_string() // '.';
    }
    else{
        say $gff;
    }
}


=head1 NAME

reverse_coords.pl - reverse coordinates of gff entries with sequence that start
with certain prefix.

=head1 SYNOPSIS

Usage examples:

 reverse_coords.pl gff-file.gff reference.fasta

Specify prefix, which defaults to 'RC_'.

 reverse_coords.pl --prefix RC_ gff-file.gff reference.fasta

Force a certain value in column 7

 reverse_coords.pl --plus  gff-file.gff reference.fasta
 reverse_coords.pl --minus gff-file.gff reference.fasta
 reverse_coords.pl --dot   gff-file.gff reference.fasta

=head1 OPTIONS

=over

=item  --prefix <prefix>

Only reverse entries when the sequence has this prefix.  Default: 'RC_'.
Case insensitive.

=item --plus

Force reverse strands to have a '+' in column 7

=item --minus

Force reverse strands to have a '-' in column 7

=item --dot

Force reverse strands to have a '.' in column 7
=back

=cut
