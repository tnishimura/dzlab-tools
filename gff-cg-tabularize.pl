#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;

use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;

use Pod::Usage;
use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "output|o=s" => \(my $output),
    "cg|c=s" => \(my $cg_file),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || ! $output || !$cg_file);  

my $p = GFF::Parser->new(file => $cg_file, normalize => 0);
open my $out, '>', $output;

while (1){
    my $gff1 = $p->next;
    my $gff2 = $p->next;
    if ((defined $gff1 && ! defined $gff2) || (!defined $gff1 && defined $gff2)){
        die "uneven # of lines?";
    }
    if (!defined $gff1 && ! defined $gff2){
        last;
    }
    my $seqid1 = $gff1->sequence;
    my $seqid2 = $gff2->sequence;
    my $start1 = $gff1->start;
    my $start2 = $gff2->start;
    my $end1 = $gff1->end;
    my $end2 = $gff2->end;
    my $c1 = $gff1->get_attribute('c') // 0;
    my $c2 = $gff2->get_attribute('c') // 0;
    my $t1 = $gff1->get_attribute('t') // 0;
    my $t2 = $gff2->get_attribute('t') // 0;

    if ($start1 != $end1 || $start2 != $end2 ){
        die "not a single-c file?";
    }

    if (lc $seqid1 ne lc $seqid2 || $start1 != $start2-1){
        die "cg file does not contain adjacent positions? see line $.";
    }

    say $out join "\t", lc $seqid1, $start1, $start2, $c1, $t1, $c2, $t2;
}
close $out;

=head1 gff-cg-tabularize.pl 

Usage examples:

 gff-cg-tabularize.pl -o output.txt input-cg-single-c.gff 

=cut

