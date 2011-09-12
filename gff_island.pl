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
use Log::Log4perl qw/:easy/;
Log::Log4perl->easy_init( { 
    level    => $DEBUG,
    #file     => ">run.log",
    layout   => '%d{HH:mm:ss} %p> (%L) %M - %m%n',
} );
my $logger = get_logger();

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || !$opt_input || !$opt_output;

$logger->info("threshold: $opt_threshold");
$logger->info("  max gap: $opt_max_gap");

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}

my $p = GFF::Parser->new(file => $opt_input);

sub land{
    my $gff = shift;
    #if (defined $gff->score){
    #die Dumper $gff;
    #}
    return defined $gff->score && $gff->score >= $opt_threshold;
}
sub unknown{
    my $gff = shift;
    return ! defined $gff->score;
}
sub sea{
    my $gff = shift;
    return defined $gff->score && $gff->score < $opt_threshold;
}

my $previous;
my ($c, $t, $n, $start, $end, $gap) = (0,0,0,undef,0,0);
while (my $gff = $p->next()){
    my $current_c;
    my $current_t;
    my $current_n;
    unless ($opt_just_scores){
        $current_c = $gff->get_attribute('c');
        $current_t = $gff->get_attribute('t');
        $current_n = $gff->get_attribute('n');
    }
    my $current_start = $gff->start();
    my $current_end = $gff->end();

    if (!$current_end || !$current_start){
        die sprintf "%s has no start or end", $gff->to_string;
    } 
    if (defined $start && $start>$current_start){
        if ($previous->sequence ne $gff->sequence){
            # okay, sequence changed.
            if ($start){
                blit($gff,$start,$end,$c,$t);
                ($c, $t, $n, $gap) = (0,0,0,0);
                undef $start;
                undef $end;
            }
        }else{
            die "gff not sorted????";
        }
    }
        
    # case 1: no on island, and we got land
    if (!$start && land($gff)){
        $start = $current_start;
        $end   = $current_end;
        unless ($opt_just_scores){
            $c += $current_c;
            $t += $current_t;
            $n += $current_n;
        }
    } 
    # case 2: not on island, no land
    elsif (!$start && ! land($gff)){
        # do nothing
    }
    # case 3: on island, land
    elsif ($start && land($gff)){
        # extend island
        $gap = 0;
        $end = $current_end;
        unless ($opt_just_scores){
            $c += $current_c;
            $t += $current_t;
            $n += $current_n;
        }
    }
    # case 4: on island, unknown
    elsif ($start && unknown($gff)){
        if ($gap >= $opt_max_gap){
            blit($gff,$start,$end,$c,$t);
            ($c, $t, $n, $gap) = (0,0,0,0);
            undef $start;
            undef $end;
        } else {
            ++$gap;
        }
    }
    elsif ($start && sea($gff)){
        blit($gff,$start,$end,$c,$t);
        ($c, $t, $n, $gap) = (0,0,0,0);
        undef $start;
        undef $end;
    }
    $previous = $gff;
}

if ($start){
    blit($previous,$start,$end,$c,$t);
    ($c, $t, $n, $gap) = (0,0,0,0);
    undef $start;
    undef $end;
}

sub blit{
    my ($gff,$start,$end,$c,$t) = @_;

    if ($opt_just_scores){
        my $score = $end - $start + 1;

        say join "\t", $gff->sequence // '.', $gff->source // '.', $gff->feature // '.',
        $start, $end, $score, $gff->strand // '.', '.', '.';
    } else {
        my $score = $c + $t > 0 ? $c / ($c + $t) : 0;

        say join "\t", $gff->sequence, $gff->source, $gff->feature,
        $start, $end, $score, $gff->strand // '.', '.', "n=$n;c=$c;t=$t";
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
    gap.type:        int, gap >= 0 

=item  -t <threshold> | --threshold <threshold>

GFF lines with scores less than this are considered empty (part of a gap). default .01.

=for Euclid
    threshold.default:     0.01

=item -s | --just-scores

Do not collect c, t, n counts from column 9 -- just go by score.

=item -h | --help

=back

=cut


