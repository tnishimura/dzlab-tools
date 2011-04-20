#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use List::Util qw/max/;
use FindBin;
use lib "$FindBin::Bin/lib";
use DZUtil qw/overlap/;
use GFF::Slurp;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || !$opt_genes_gff || !$opt_exons_gff || !$opt_output;

if (!$opt_five_prime && !$opt_three_prime){
    $opt_five_prime = 1;
}

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}

my $genes = gff_slurp_index($opt_genes_gff,$opt_gene_locus);
my $exons = gff_slurp_index($opt_exons_gff,'*' . $opt_exon_locus);

open my $in, '<', $opt_input;

TOP:
while (defined(my $line = <$in>)){
    chomp $line;
    my ($locus, @scores) = split /\t/, $line;
    my $genes = $genes->{$locus};
    my $exons = $exons->{$locus};
    if (! defined $genes){
        warn "$locus not found in genes GFF file $opt_genes_gff? skipping...";
        next TOP;
    }
    if (! defined $exons){
        warn "$locus not found in exons GFF file $opt_exons_gff? skipping...";
        next TOP;
    }

    print $locus;

    my $gene_start = $genes->[0]->start();
    my $gene_end   = $genes->[0]->end();

    my $reverse
        = ( $opt_five_prime and $genes->[0]->strand() eq q{-} )
            || ( $opt_three_prime and $genes->[0]->strand() eq q{+} )
                || 0;

    my $current_start = $reverse ? $gene_end + $opt_distance           : $gene_start - $opt_distance;
    my $current_end   = $reverse ? $current_start - $opt_bin_width + 1 : $current_start + $opt_bin_width - 1;

    for my $score (@scores) {
        if ($score eq 'na'){
            print "\tna";
        }
        elsif (grep { 
                overlap([$current_start,$current_end], [$_->start,$_->end]) / $opt_bin_width >= ($opt_threshold / 100)
            } @$exons){
            print "\t$score";
        } 
        else {
            print "\tna";
            print "*" if $opt_debug;
        }
        print "($current_start,$current_end)" if $opt_debug;
    } continue {
        $current_start += ($reverse ? -1 : 1) * $opt_bin_width;
        $current_end   += ($reverse ? -1 : 1) * $opt_bin_width;
    }
    print "\n";
}

=head1 NAME
 
ends_filter.pl - given an ends file, a gene annotation file and an exon annotation file, remove all scores which are not in exons
from the ends file.
 
=head1 SYNOPSIS

Usage: 

Bin-width 100, centered on 5' end, 5000 distance:  

 ends_filter.pl -5 -b 100 -d 5000 -g TAIR8_genes.gff -e TAIR8_exons.gff -gl ID -el Parent -o output.ends input.ends

=head1 OPTIONS

=over

=item  -b <width> | --bin-width <width>

Bin width.  Should be the same as passed to ends_analysis.pl. Default to 100.

=for Euclid
    width.default:     100

=item  -d <distance> | --distance <distance>

Distance from end to analyze.  Should be same as passed to ends_analysis.pl. Default to 5000.

=for Euclid
    distance.default:     5000

=item -5 | --five-prime

Center analysis on 5' end. Default. Should be same as passed to ends_analysis.pl

=item -3 | --three-prime

Center analysis on 3' end. 

=item  -g <gff> | --genes-gff <gff>

Gene annotation file, should be the same as for ends_analysis.pl.

=for Euclid
    gff.type:        readable

=item  -e <gff> | --exons-gff <gff>

Exon annotation file.

=for Euclid
    gff.type:        readable

=item  -gl <locus> | --gene-locus <locus>

Locus tag in --genes-gff.  Default to 'ID', which means the script expects
something like 'ID=AT1G01010' in the attributes column (col 9).

=for Euclid
    locus.default:     'ID'

=item  -el <locus> | --exon-locus <locus>

Locus tag in --genes-gff.  Default to 'Parent', which means the script expects
something like 'ID=AT1G01010.1' in the attributes column (col 9).

=for Euclid
    locus.default:     'Parent'

=item  -t <threshold> | --threshold <threshold>

If over threshold% of a bin overlaps with an exon, keep its score. Otherwise,
replace with 'na'. Default 50, meaning if 50% of a bin overlaps with any exon,
its score is kept. 

=for Euclid
    threshold.default:     50
    threshold.type:        int, threshold >= 0 && threshold <= 100

=item  <input>

Input file.  Should be the output file of ends_analysis.pl.

=item  --debug

Print each score with position, and eliminated scores with '*'.  

=item -o <file> | --output <file>

=item -h | --help

=back

=cut

