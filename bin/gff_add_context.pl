#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
END {close STDOUT}
$| = 1;
use FindBin;
use lib "$FindBin::Bin/../lib";
use FastaReader;
use GFF::Parser;

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "reference|r=s" => \(my $reference),
    "debug|d" => \(my $debug),
);
pod2usage(-verbose => 2, -noperldoc => 1) 
if (!$result || ! $reference || ! -f $reference);

my $fr = FastaReader->new(file => $reference, slurp => 1);
my $p = GFF::Parser->new(file => \*ARGV, normalize => 0);

while (defined(my $gff = $p->next())){
    my $seq = $gff->sequence();
    my $start = $gff->start();
    my $end = $gff->end();
    die "$0 only supports single-c files" if $start != $end;
    my $context = $fr->get_context($seq, $start, undef_on_nonc => 1);
    if ($context){
        $gff->feature($context);
    }
    say $gff;
}


=head1 NAME

gff_add_context.pl - add methylation context for all C/G positions in GFF file.


=head1 SYNOPSIS

Only supports single-c files.  Non C/G positions are untouched.

 gff_add_context.pl -r genome.fasta single-c.gff

=cut

