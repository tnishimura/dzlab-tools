#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use GFF;

use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) if !$opt_output;

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh;
}
my $p = GFF::Parser->new(file => $opt_infile);

while (defined(my $gff = $p->next())){
    say join "\t", 
    $gff->sequence // '.', $gff->source // '.', $gff->feature // '.', $gff->start // '.', $gff->end // '.', 
    $gff->score // '.', 
    $gff->strand // '.', 
    ($gff->get_column('c')//0) + ($gff->get_column('t')//0), 
    $gff->attribute_string // '.';
}

=head1 NAME

fracmeth_sum.pl - Sum the 'c' and 't' in column 9 and put it in column 8.

=head1 SYNOPSIS

Usage examples:

 fracmeth_sum.pl -o output.gff file.gff

 fracmeth_sum.pl -o - file.gff

=head1 OPTIONS

=over

=item --help | -h

=item  -o <file> | --output <file>

=item <infile>

=for Euclid
    infile.type:        readable

=back

=cut

