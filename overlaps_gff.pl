#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use lib "$FindBin::Bin/Tree-Range/lib";
use GFF::Tree;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;


pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || !$opt_gff || !$opt_input;

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}

my $gt = GFF::Tree->new(file => $opt_gff);

my $p = GFF::Parser->new(file => $opt_input);

while (defined(my $gff = $p->next())){
    my ($seq, $start, $end) = ($gff->sequence, $gff->start, $gff->end);
    my @results = $gt->search_overlap($seq,$start,$end);
    my $outputted;

    for my $result (@results){
        my $overlap = $result->{overlap};
        if ($opt_proportion_threshold == 0 || 
            ($opt_proportion_threshold == 1 && $overlap == $gff->length) || 
            $overlap / $gff->length() >= $opt_proportion_threshold
        ){
            if (defined $opt_tag && defined(my $locus = $result->{item}->get_column($opt_tag))){
                if ($opt_tag_only){
                    $gff->attribute_string($locus);
                    say $gff;
                }
                else{
                    say "$gff;$opt_tag=$locus";
                }
                $outputted = 1;
            }
            else{
                say STDERR "Warning: overlapping but no $opt_tag for " . $result->{item}; 
            }
        }
    }

    if (! $outputted && $opt_no_skip){
        $gff->score(undef);
        $gff->attribute_string(undef) if $opt_tag_only;
        say $gff;
    }
}


=head1 NAME

overlaps_gff.pl 

=head1 SYNOPSIS

Usage examples:

 overlaps_gff.pl -g annotation.gff -k -t 'ID' -o output.gff input.gff

=head1 REQUIRED ARGUMENTS

=over

=back

=head1 OPTIONS

=over


=item  -g <file> | --gff <file>

Annotation file.

=for Euclid
    file.type:        readable

=item  <input>

=for Euclid
    input.type:        readable

=item  -k | --no-skip

Print line from input.gff even if no overlap is found.

=item  -t <tag> | --tag <tag>

Locus tag in --gff annotation file. Optional. Default to giving you the entire attribute string (col 9) from the annotation.

=item -o <file> | --output <file>

=for Euclid
    file.default:     '-'

=item  -p <thresh> | --proportion-threshold <thresh>

Input GFF lines needs at least this proportion covered by annotation GFF. For example,
if this is set to .50, and the input GFF is windowed by 50, at least 25 bases need to be overlapping
with an annotation in order to be reported.  If not, skipped (unless --no-skip is used, in which case it's 
reported with a null score). Default to 0, mean any overlap is enough.

=for Euclid
    thresh.default:     0
    thresh.type:        number, thresh >= 0 && thresh <= 1
    thresh.type.error:  <thresh> must be between 0 and 1

=item  -to  | --tag-only 

Output only tag in column 9, instead of original col 9 + tag.

=item --help | -h

=back

=cut



