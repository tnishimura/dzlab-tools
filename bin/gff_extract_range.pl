#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/lib";
use DZUtil qw/overlap/;
use GFF::Parser;
use FastaReader;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if ($opt_help || ! (scalar keys %opt_range) || ! ($opt_range{start} <= $opt_range{end} ));

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}
my $p = GFF::Parser->new(file => $opt_input);
my $reverse_coord;
if ($opt_reverse_coord){
    $reverse_coord = FastaReader->new(file => $opt_reverse_coord, slurp => 0);
}


LOOP:
while (defined(my $gff = $p->next())){
    if (defined $opt_sequence){
        next LOOP if ! defined $gff->sequence();
        next LOOP if uc($gff->sequence()) ne uc($opt_sequence);
    }
    if (defined $gff->strand() and $gff->strand() eq '-' and defined $reverse_coord){
        my $end   = $reverse_coord->reverse2forward($gff->sequence(),$gff->start());
        my $start = $reverse_coord->reverse2forward($gff->sequence(),$gff->end());
        $gff->start($start);
        $gff->end($end);
    }

    if ( $gff->start <= $gff->end ){
        #die "asdf" if $counter++ == 100; 
        my $overlap   = overlap([$gff->start, $gff->end], [$opt_range{start}, $opt_range{end}]);
        my $targetlen = $gff->end - $gff->start + 1;
        my $querylen  = $opt_range{end} - $opt_range{start} + 1;

        
        if (($opt_threshold == 0   && $overlap) || # handle extreme cases separately so no rounding errors
            ($opt_threshold == 100 && $overlap == $targetlen) || 
            ($overlap / $targetlen > ($opt_threshold / 100))){
            if  (!$opt_invert){
                say $gff->to_string();
            }
        }
        elsif ($opt_invert){
            say $gff->to_string();
        }
    }

}

=head1 NAME
 
gff_extract_range.pl
 
=head1 SYNOPSIS

 gff_extract_range.pl

=head1 OPTIONS

=over

=item  -s <seqid> | --sequence <seqid>

If given, match sequence ID (column 1) as well as range. If not given
(default), return all matching ranged regardless of sequence.  If you are
working with a single-chromosome file, you don't need to specify this.

=item  -r <start> <end> | --range <start> <end> 

Query range.

=for Euclid
    start.type:      int, start >= 1 
    end.type:        int, end >= 1 

=item  -t <percentage> | --threshold <percentage>

Percent of target window which needs to be covered by the query range in order
to be reported.  If set to 0, any target windows with any overlap will be
reported. Default 100.

=for Euclid
    percentage.default: 100
    percentage.type:    int, percentage >= 0 && percentage <= 100

=item -o <file> | --output <file>

=for Euclid
    file.default: '-'

=item  -i  | --invert

Print lines that do NOT match the filtering criteria.

=item  -rc <reference_genome> | --reverse-coord <reference_genome>

Fix the the case where minus strand coordinates are with respect to the 3'
end.  This option usually only makes sense for .eland3.post files, from the
bs-seq(uel) pipelines.  You have to pass a reference genome b/c it needs to
know the length of chromosomes. Perhaps this should be separate script...

=for Euclid
    reference_genome.type:        readable

=item <input>

=for Euclid
    input.type: readable

=item -h | --help

=back

=cut

