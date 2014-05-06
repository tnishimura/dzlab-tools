#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/lib";
use FastaReader;

use Pod::Usage;
use Getopt::Long;

my $reference = shift or pod2usage(-verbose => 2, -noperldoc => 1);
my $output_file_cg  = "$reference.cg-sites.txt";
my $output_file_chg = "$reference.chg-sites.txt";
my $output_file_chh = "$reference.chh-sites.txt";

open my $fh_cg, '>', $output_file_cg;
open my $fh_chg, '>', $output_file_chg;
open my $fh_chh, '>', $output_file_chh;

my $fr = FastaReader->new(file => $reference, slurp => 1);

my $counter = 0;
for my $seqid ($fr->sequence_list) {
    warn "indexing $seqid";
    my $seqid = lc($seqid);
    my $seq = $fr->get($seqid);
    while ($seq =~ m{CG}g){
        my $pos1 = $-[0] + 1;
        my $pos2 = $-[0] + 2;
        say $fh_cg join "\t", $seqid, $pos1, "+", "CG";
        say $fh_cg join "\t", $seqid, $pos2, "-", "CG";
    }

    while ($seq =~ m{(?=C[ACT]G)}g){
        my $pos1 = $-[0] + 1;
        say $fh_chg join "\t", $seqid, $pos1, "+", "CHG";
    }

    while ($seq =~ m{(?=C[ATG]G)}g){
        my $pos2 = $-[0] + 3;
        say $fh_chg join "\t", $seqid, $pos2, "-", "CHG";
    }

    while ($seq =~ m{(?=C[ACT][ACT])}g){
        my $pos1 = $-[0] + 1;
        say $fh_chh join "\t", $seqid, $pos1, "+", "CHH";
    }

    while ($seq =~ m{(?=[ATG][ATG]G)}g){
        my $pos2 = $-[0] + 3;
        say $fh_chh join "\t", $seqid, $pos2, "-", "CHH";
    }
}
close $fh_cg;
close $fh_chg;
close $fh_chh;

=head1 methyl-sites.pl 

 methyl-sites.pl referencee.fas

=cut

