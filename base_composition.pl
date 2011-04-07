#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF;
use GFF::Parser;
use Fasta;

use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if ($opt_help || !$opt_gff || !$opt_input || !$opt_output);

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}

my $fasta = slurp_fasta($opt_input);
my $parser = GFF::Parser->new(skip => 1,file => $opt_gff);

while (my $gff = $parser->next()){
    my $body = fasta_subseq($fasta,$gff->sequence,$gff->start, $gff->end, coord => 'f', rc => $gff->strand eq '-' ? 1 : 0);
    my $len = $gff->end - $gff->start +1;
    my $a = ($body=~tr/aA/aA/) / $len;
    my $c = ($body=~tr/cC/cC/) / $len;
    my $g = ($body=~tr/gG/gG/) / $len;
    my $t = ($body=~tr/tT/tT/) / $len;
    say join "\t",
    $gff->sequence, 
    $gff->source,
    $gff->feature,
    $gff->start,
    $gff->end,
    $gff->score // 0, 
    $gff->strand // '+',
    $gff->frame // 0,
    sprintf("a=%.5f;c=%.5f;g=%.5f;t=%.5f",$a,$c,$g,$t) 
    . (defined $gff->attribute_string ? ";" . $gff->attribute_string : "");
}
if ($opt_output ne '-'){
    close \*STDOUT;
}

=head1 NAME

base_composition.pl - Given a gff annotation file and a corresponding reference fasta file, add the approrpiate ratios
for each base to the attribute field (column 9) of the gff and spit it back out.

=head1 SYNOPSIS

Usage examples:

 base_composition.pl -g annotation.gff -o annotation-with-count.gff reference.fasta 

=head1 OPTIONS

=over

=item  -g <file> | --gff <file>

GFF annotation file.  

=for Euclid
    file.type:        readable

=item  <input> 

Input fasta file.

=for Euclid
    input.type:        readable

=item  -o <file> | --output <file>

output. Use '-' (without quotes) to print to screen.

=item --help | -h

=back

=cut




