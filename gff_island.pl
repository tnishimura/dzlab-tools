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
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || !$opt_input || !$opt_output;

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}

my $p = GFF::Parser->new(file => $opt_input);

sub land{
    my $gff = shift;
    return defined $gff->score && $gff->score >= $opt_threshold;
}

my ($c, $t, $n, $start, $gap) = (0,0,0,undef,0);
while (my $gff = $p->next()){
    my $current_c = $gff->get_attribute('c');
    my $current_t = $gff->get_attribute('t');
    my $current_n = $gff->get_attribute('n');
    my $current_start = $gff->start();
    my $current_end = $gff->end();

    if (!$current_end || !$current_start){
        die sprintf "%s has no start or end", $gff->to_string;
    } 
    if (defined $start && $start>$current_start){
        die "gff not sorted????";
    }
        

    # case 1: no on island, and we got land
    if (!$start && land($gff)){
        $start = $current_start;
        $c += $current_c;
        $t += $current_t;
        $n += $current_n;
    } 
    # case 2: not on island, no land
    elsif (!$start && ! land($gff)){
        # do nothing
    }
    # case 3: on island, land
    elsif ($start && land($gff)){
        # extend island
        $c += $current_c;
        $t += $current_t;
        $n += $current_n;
    }

    elsif ($start && ! land($gff)){
        if ($gap > $opt_max_gap){
            my $score = $c + $t > 0 ? $c / ($c + $t) : 0;
            say join "\t", 
            $gff->sequence, 
            $gff->source,
            $gff->feature,
            $start,
            $current_end,
            $score,
            $gff->strand // '.',
            '.',
            "n=$n;c=$c;t=$t";
            ($c, $t, $n, $gap) = (0,0,0,0);
            undef $start;
        } else {
            ++$gap;
        }
    }
}
if ($opt_output ne '-'){
    close \*STDOUT;
}

=head1 NAME
 
gff_island.pl - find islands of 
 
=head1 SYNOPSIS

 gff_island.pl

=head1 OPTIONS

=over

=item  -i <file> | --input <file>

=for Euclid
    file.type:        readable

=item -o <file> | --output <file>

=item  -g <gap> | --max-gap <gap>

Maximum number of GFF lines (does NOT take into account window width) that a single island is allowed to extend to.
Default 100.

=for Euclid
    gap.default:     100
    gap.type:        int, gap >= 1 

=item  -t <threshold> | --threshold <threshold>

GFF lines with scores less than this are considered empty (part of a gap). default .01.

=for Euclid
    threshold.default:     0.01


=item -h | --help

=back

=cut


