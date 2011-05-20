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

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if ($opt_help || ! (scalar keys %opt_range) || ! ($opt_range{start} <= $opt_range{end} ));

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}
my $p = GFF::Parser->new(file => $opt_input);

#my $counter = 0;
while (defined(my $gff = $p->next())){
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

=item <input>

=for Euclid
    input.type: readable

=item -h | --help

=back

=cut


