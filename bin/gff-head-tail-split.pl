#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;

use FindBin;
use lib "$FindBin::Bin/../lib";
use GFF::Parser;

use Pod::Usage;
use Getopt::Long;
use List::Util qw/max min/;


my $result = GetOptions (
    "file|f=s" => \(my $file),
    "size|s=i" => \(my $size = 500),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || ! $file);  

my $head_file = $file =~ s/\.gff$/-head-$size.gff/ir;
my $tail_file = $file =~ s/\.gff$/-tail-$size.gff/ir;
my $head = IO::File->new($head_file, 'w');
my $tail = IO::File->new($tail_file, 'w');

my $p = GFF::Parser->new(file => $file, normalize => 0);
while (defined(my $gff = $p->next())){
    my $start = $gff->start;
    my $end   = $gff->end;
    my $feature = $gff->feature;

    # forward
    if (! $gff->is_reverse){
        my $mid = min($start + $size - 1, $end);
        $gff->start($start);
        $gff->end($mid);
        $gff->feature("$feature-head");
        $head->print($gff . "\n");

        if ($mid < $end){
            $gff->start($mid + 1);
            $gff->end($end);
            $gff->feature("$feature-tail");
            $tail->print($gff . "\n");
        }
    }
    # reverse
    else{
        my $mid = max($end - $size + 1, $start);
        $gff->start($mid);
        $gff->end($end);
        $gff->feature("$feature-head");
        $head->print($gff . "\n");

        if ($mid > $start){
            $gff->start($start);
            $gff->end($mid - 1);
            $gff->feature("$feature-tail");
            $tail->print($gff . "\n");
        }
    }
}
=head1 gff-head-tail-split.pl 

 gff-head-tail-split.pl -s 500 -f annotation.gff

Split annotation entries into a 'head' of -s bp (default: 500) and a tail of (length - size) bp.

(it turns this:

 chr1	TAIR10	gene	3631	5899	.	+	.	ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
 chr1	TAIR10	gene	5928	8737	.	-	.	ID=AT1G01020;Note=protein_coding_gene;Name=AT1G01020
 chr1	TAIR10	gene	11649	13714	.	-	.	ID=AT1G01030;Note=protein_coding_gene;Name=AT1G01030
 chr1	TAIR10	gene	23146	31227	.	+	.	ID=AT1G01040;Note=protein_coding_gene;Name=AT1G01040
 chr1	TAIR10	gene	28500	28706	.	+	.	ID=AT1G01046;Note=miRNA;Name=AT1G01046
 chr1	TAIR10	gene	31170	33153	.	-	.	ID=AT1G01050;Note=protein_coding_gene;Name=AT1G01050
 chr1	TAIR10	gene	33379	37871	.	-	.	ID=AT1G01060;Note=protein_coding_gene;Name=AT1G01060
 chr1	TAIR10	gene	38752	40944	.	-	.	ID=AT1G01070;Note=protein_coding_gene;Name=AT1G01070
 chr1	TAIR10	gene	44677	44787	.	+	.	ID=AT1G01073;Note=protein_coding_gene;Name=AT1G01073
 chr1	TAIR10	gene	45296	47019	.	-	.	ID=AT1G01080;Note=protein_coding_gene;Name=AT1G01080

into this:

 chr1	TAIR10	gene-head	3631	4130	.	+	.	ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
 chr1	TAIR10	gene-head	8238	8737	.	-	.	ID=AT1G01020;Note=protein_coding_gene;Name=AT1G01020
 chr1	TAIR10	gene-head	13215	13714	.	-	.	ID=AT1G01030;Note=protein_coding_gene;Name=AT1G01030
 chr1	TAIR10	gene-head	23146	23645	.	+	.	ID=AT1G01040;Note=protein_coding_gene;Name=AT1G01040
 chr1	TAIR10	gene-head	28500	28706	.	+	.	ID=AT1G01046;Note=miRNA;Name=AT1G01046
 chr1	TAIR10	gene-head	32654	33153	.	-	.	ID=AT1G01050;Note=protein_coding_gene;Name=AT1G01050
 chr1	TAIR10	gene-head	37372	37871	.	-	.	ID=AT1G01060;Note=protein_coding_gene;Name=AT1G01060
 chr1	TAIR10	gene-head	40445	40944	.	-	.	ID=AT1G01070;Note=protein_coding_gene;Name=AT1G01070
 chr1	TAIR10	gene-head	44677	44787	.	+	.	ID=AT1G01073;Note=protein_coding_gene;Name=AT1G01073
 chr1	TAIR10	gene-head	46520	47019	.	-	.	ID=AT1G01080;Note=protein_coding_gene;Name=AT1G01080

and this:

 chr1	TAIR10	gene-tail	4131	5899	.	+	.	ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
 chr1	TAIR10	gene-tail	5928	8237	.	-	.	ID=AT1G01020;Note=protein_coding_gene;Name=AT1G01020
 chr1	TAIR10	gene-tail	11649	13214	.	-	.	ID=AT1G01030;Note=protein_coding_gene;Name=AT1G01030
 chr1	TAIR10	gene-tail	23646	31227	.	+	.	ID=AT1G01040;Note=protein_coding_gene;Name=AT1G01040
 chr1	TAIR10	gene-tail	31170	32653	.	-	.	ID=AT1G01050;Note=protein_coding_gene;Name=AT1G01050
 chr1	TAIR10	gene-tail	33379	37371	.	-	.	ID=AT1G01060;Note=protein_coding_gene;Name=AT1G01060
 chr1	TAIR10	gene-tail	38752	40444	.	-	.	ID=AT1G01070;Note=protein_coding_gene;Name=AT1G01070
 chr1	TAIR10	gene-tail	45296	46519	.	-	.	ID=AT1G01080;Note=protein_coding_gene;Name=AT1G01080
 chr1	TAIR10	gene-tail	47485	48786	.	-	.	ID=AT1G01090;Note=protein_coding_gene;Name=AT1G01090
 chr1	TAIR10	gene-tail	50075	50699	.	-	.	ID=AT1G01100;Note=protein_coding_gene;Name=AT1G01100

)

=cut

