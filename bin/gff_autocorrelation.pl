#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use List::MoreUtils qw/all/;
use File::Basename qw/basename/;
use PDL::Lite;
use PDL::Stats;
use Scalar::Util qw/looks_like_number/;

use FindBin;
use lib "$FindBin::Bin/../lib";
use Launch;
use GFF::Parser;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || !$opt_gff || !$opt_input;

### intermediate file
my $opt_overlaps_file = basename($opt_input) . ".overlaps-vs-" . basename($opt_gff, '.gff');
my $sorted = "$opt_overlaps_file.sorted";


say STDERR "\$opt_overlaps_file = $opt_overlaps_file";
launch("overlaps_gff.pl -t $opt_tag -to -g $opt_gff -o ?? $opt_input", expected => $opt_overlaps_file);
launch("sort -k9,9 $opt_overlaps_file -o ??", expected => $sorted);

say STDERR "making position map";
my $posmap = make_position_map($opt_gff);
say STDERR "done making position map";

exit 0 if $opt_overlaps_only;

#######################################################################
# output
$opt_output //= $opt_overlaps_file . ".acf";
open my $fh, '>', $opt_output;

say STDERR "\$opt_output        = $opt_output";

my $parser = GFF::Parser->new(file => $sorted);

my $last_tag;
my $last_template;

while (defined(my $gff = $parser->next())){
    my ($seq, $score, $start, $end, $tag) = ($gff->sequence, $gff->score(), $gff->start, $gff->end, $gff->attribute_string);
    $score//=0;
    #say STDERR "$seq, $score, $start, $end, $tag";
    die "not single-nuc" if ($start != $end);
    die "no tag" if ! defined $tag;

    $last_tag //= $tag;
    $last_template //= get_template($posmap, $tag);

    if ($last_tag ne $tag){
        if (defined(my $auto = autocorrelation($last_template))){
            format_line($last_tag, @$auto);
        }

        $last_tag = $tag;
        $last_template = get_template($posmap, $tag);
    }

    $last_template->[get_relative_position($posmap, $tag, $start)] += $score;
}

if (defined(my $auto = autocorrelation($last_template))){
    format_line($last_tag, @$auto);
}

close $fh;

#######################################################################
# average

{
    open my $input_fh, '<', $opt_output;

    my @indices = (0 .. $opt_cutoff-1);
    my @cumulative_values = (0) x $opt_cutoff;
    my @counts = (0) x $opt_cutoff;

    while (defined(my $line = <$input_fh>)){
        chomp $line;
        my (undef, @split) = split /\t/, $line;
        @split = @split[0 .. $opt_cutoff - 1];

        for my $i (@indices){
            my $v = $split[$i];
            if (defined $v && $v ne 'nan' && looks_like_number $v){
                $cumulative_values[$i]+=$v;
                $counts[$i]++;
            }
        }
    }
    close $input_fh;

    open my $average_out_fh, '>', "$opt_output.avg";
    for my $i (@indices) {
        my $c = $counts[$i];
        my $v = $cumulative_values[$i];

        say $average_out_fh join "\t", $i, $c ? $v / $c : 0, $v, $c;
    }
    close $average_out_fh;
}

#######################################################################

sub format_line{
    my ($tag, @vals) = @_;
    say $fh join "\t", $tag, map {defined($_) ? sprintf("%.7f", $_) : '0'} @vals[0 .. $opt_cutoff - 1];
}

sub make_position_map{
    my $file = shift;
    my $p = GFF::Parser->new(file => $file);
    my %accum;
    my $counter = 0;
    while (defined(my $gff = $p->next)){
        if (defined(my $t = $gff->get_column($opt_tag))){
            $accum{$t} = [$gff->strand(), $gff->start(), $gff->end()];
        }
        say STDERR $counter if ++$counter % 1000 == 0;
    }
    return \%accum;
}

sub get_template{
    my ($posmap, $tag) = @_;
    my ($strand, $start, $end) = @{$posmap->{$tag}};
    return [(0) x ($end - $start + 1)];
}

sub get_relative_position{
    my ($posmap, $tag, $position) = @_;
    my ($strand, $start, $end) = @{$posmap->{$tag}};
    if ($start <= $position && $position <= $end){
        if ($strand eq '+'){
            return $position - $start;
        }
        elsif ($strand eq '-'){
            return $end - $position;
        }
        else {
            die "impossible 1";
        }
    }
    else{
        die "impossible 2";
    }
}

sub autocorrelation{
    my $points = shift;
    if (@$points == 0 || all { $_ == $points->[0] } @$points){
        return;
    }
    if (@$points > $opt_cutoff){
        @$points = @{$points}[0 .. $opt_cutoff - 1];
    }
    my $pdl = pdl $points;

    return [$pdl->acf()->list];
}


=head1 NAME

overlaps_gff.pl 

=head1 SYNOPSIS

Usage examples:

 overlaps_gff.pl -g annotation.gff -k -t 'ID' -o output.gff input.gff

=head1 OPTIONS

=over

=item  -g <file> | --gff <file>

Annotation file.

=for Euclid
    file.type:        readable

=item  <input>

=for Euclid
    input.type:        readable

=item  -t <tag> | --tag <tag>

Locus tag in --gff annotation file. Optional. Default to giving you the entire
attribute string (col 9) from the annotation.

=for Euclid
    tag.default:     'ID'

=item  -o <file> | --output <file>

=item --help | -h

=item  -c <length> | --cutoff <length>

=for Euclid
    length.default:     2000
    length.type:        int, length >= 1 && length <= 10000

=item  -oo | --overlaps-only 

=back

=cut

