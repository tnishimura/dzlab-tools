#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use Scalar::Util qw/looks_like_number/;


END {close STDOUT}
$| = 1;


use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) if $opt_help;

my $p = GFF::Parser->new(file => $opt_input);
open my $outfh, '>', $opt_output;

while (defined(my $gff = $p->next())){
    my @gff_split = (
        $gff->sequence // '.', 
        $gff->source // '.', 
        $gff->feature // '.', 
        undef,
        undef,
        $gff->score // 1, 
        $gff->strand // '.', 
        $gff->frame // '.', 
        $gff->attribute_string // '.'); 
    if (looks_like_number $gff->start && looks_like_number $gff->end && $gff->start < $gff->end){
        for my $index ($gff->start .. $gff->end){
            @gff_split[3,4] = ($index, $index);
            say $outfh join "\t", @gff_split;
        }
    }
}



=head1 NAME

gff_explode.pl - explode windows into single-nucleotide positions

=head1 SYNOPSIS

Usage examples:

 gff_explode.pl [options]...

=head1 REQUIRED ARGUMENTS

=over

=item  -i <file> | --input <file>

=for Euclid
    file.type:        readable

=item  -o <file> | --output <file>

=back

=head1 OPTIONS

=over

=item --help | -h

=back

=cut
