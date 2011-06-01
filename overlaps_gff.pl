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
    my $gffstring = $gff->to_string;

    if (@results){
        for my $result (@results){
            if (my $locus = $result->{item}->get_column($opt_tag)){
                say "$gffstring;ID=$locus";
            }
            else{
                say "$gffstring;" . $result->{item}->attribute_string;
            }
        }
    }
    elsif ($opt_no_skip){
        $gff->score(undef);
        say $gff->to_string;
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

=for Euclid
    file.type:        readable

=item  <input>

=for Euclid
    input.type:        readable

=item  -k | --no-skip

Print line from input.gff even if no overlap is found.

=item  -t <tag> | --tag <tag>

Locus tag in --gff annotation file. Optional, default to 'ID'. 

=for Euclid
    tag.default:     'ID'

=item -o <file> | --output <file>

=for Euclid
    file.default:     '-'

=item --help | -h

=back

=cut



