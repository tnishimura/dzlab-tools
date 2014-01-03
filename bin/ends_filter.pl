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
my $counter = 0;

TOP:
while (defined(my $line = <$in>)){ # for each line in ends
    chomp $line;
    my ($locus, @scores) = split /\t/, $line;
    my $genes = $genes->{$locus};
    my $exons = $exons->{$locus};
    #say STDERR join "\n", @$exons;
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

    my $strand = $genes->[0]->strand() eq q{+} ? 'f' : $genes->[0]->strand() eq q{-} ? 'r' : die "wha?";
    my $prime = $opt_five_prime ? 5 : $opt_three_prime ? 3 : die "wha?";

    #say STDERR "$opt_input $strand $prime";

    my $current_start;
    my $current_end;
    if ($prime == 5 && $strand eq 'f'){
        $current_start = $gene_start - $opt_distance;
        $current_end   = $current_start + $opt_bin_width - 1;
    }
    elsif ($prime == 3 && $strand eq 'f'){
        $current_start = $gene_end - $opt_distance;
        $current_end   = $current_start + $opt_bin_width - 1;
    }
    elsif ($prime == 5 && $strand eq 'r'){
        $current_start = $gene_end + $opt_distance;
        $current_end   = $current_start - $opt_bin_width + 1;
    }
    elsif ($prime == 3 && $strand eq 'r'){
        $current_start = $gene_start + $opt_distance;
        $current_end   = $current_start - $opt_bin_width + 1;
    }

    my $e;
    #my $xcounter = 0;
    for my $score (@scores) {
        $e = 0;
        if ($score eq 'na'){
            print "\tna";
        }
        # If it overlaps with exon, always allow it.
        elsif (grep { 
                overlap([sort {$a<=>$b} $current_start,$current_end], [$_->start,$_->end]) / $opt_bin_width >= ($opt_threshold / 100)
            } @$exons){
            print "\t$score";
            $counter++;
            $e = 1;
        } 
        # otherwise, if we're not filtering outside and current is not in gene, allow it
        elsif (!$opt_everywhere && overlap([sort {$a<=>$b} $current_start, $current_end],[$gene_start,$gene_end]) == 0){
            print "\t$score";
            # BUG: if a bin doesn't overlap sufficiently with an exon (previous block)
            # and does overlap with gene even 1 bp, it'll be eliminated...
        }
        # else eliminate
        else {
            print "\tna";
            print "*" if $opt_debug;
        }
        print "($current_start,$current_end)" if $opt_debug;
        #say STDERR "$current_start -> $current_end ($e) #" . $xcounter++;
        $current_start += ($strand eq 'r' ? -1 : 1) * $opt_bin_width;
        $current_end   += ($strand eq 'r' ? -1 : 1) * $opt_bin_width;
    }
    print "\n";
}

#say STDERR $counter;

# Chr1    TAIR8   gene    3631    5899    .       +       .       ID=AT1G01010
# Chr1    TAIR8   exon    3631    3913    .       +       .       Parent=AT1G01010.1
# Chr1    TAIR8   exon    3996    4276    .       +       .       Parent=AT1G01010.1
# Chr1    TAIR8   exon    4486    4605    .       +       .       Parent=AT1G01010.1
# Chr1    TAIR8   exon    4706    5095    .       +       .       Parent=AT1G01010.1
# Chr1    TAIR8   exon    5174    5326    .       +       .       Parent=AT1G01010.1
# Chr1    TAIR8   exon    5439    5899    .       +       .       Parent=AT1G01010.1
# Chr1    TAIR8   gene    6790    8737    .       -       .       ID=AT1G01020
# Chr1    TAIR8   exon    6790    7069    .       -       .       Parent=AT1G01020.1
# Chr1    TAIR8   exon    7157    7232    .       -       .       Parent=AT1G01020.1
# Chr1    TAIR8   exon    7384    7450    .       -       .       Parent=AT1G01020.1
# Chr1    TAIR8   exon    7564    7649    .       -       .       Parent=AT1G01020.1
# Chr1    TAIR8   exon    7762    7835    .       -       .       Parent=AT1G01020.1
# Chr1    TAIR8   exon    7942    7987    .       -       .       Parent=AT1G01020.1
# Chr1    TAIR8   exon    8236    8325    .       -       .       Parent=AT1G01020.1
# Chr1    TAIR8   exon    8417    8464    .       -       .       Parent=AT1G01020.1
# Chr1    TAIR8   exon    8571    8737    .       -       .       Parent=AT1G01020.1


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

=item  -n  | --everywhere 

By default, this script will not affect bins that are outside genes, which are clearly outside exons.  
Using this options removes scores outside genes.

=item  --debug

Print each score with position, and eliminated scores with '*'.  

=item -o <file> | --output <file>

=item -h | --help

=back

=cut

