#!/usr/bin/env perl
# oneshot for Yvonne
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


my $result = GetOptions (
    "annotation|a=s" => \(my $annotation),
    "locus-tag|t=s" => \(my $locus_tag = 'ID'),
);
if (! $result || ! $annotation || (! @ARGV && -t STDIN)){
    pod2usage(-verbose => 2, -noperldoc => 1);
}

my $annotation_parser = GFF::Parser->new(file => $annotation);
my $input_parser = GFF::Parser->new(file => \*ARGV);

my %strand;
while (defined(my $anno = $annotation_parser->next())){
    $strand{$anno->get_column($locus_tag)} = $anno->strand();
}

while (defined(my $in = $input_parser->next())){
    my $locus = $strand{$in->get_column($locus_tag)};
    if ($locus){
        $in->strand($locus);
    }
    say $in;
}

=head1 gff_add_strand_from.pl 

Usage examples:

 gff_add_strand_from.pl -a TAIR8_genes.gff input.gff > output.gff
 gff_add_strand_from.pl -a Rice_genes.gff -t Alias input.gff > output.gff

=cut
