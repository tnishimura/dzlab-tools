#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;

use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use lib "$FindBin::Bin/Tree-Range/lib";
use GFF::Tree;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use Log::Log4perl qw/:easy/;
use Counter;

Log::Log4perl->easy_init({ level => $DEBUG, layout => '%d{HH:mm:ss} %.1p > %m%n' });
my $logger = get_logger();

pod2usage(-verbose => 2,-noperldoc => 1)
if $opt_help || !$opt_gff || !$opt_input || !$opt_tag;

#######################################################################

my %do_weighting;
if (defined $opt_weigh){
    if ($opt_weigh =~ /[^ctn6]/){
        die "--weigh / -w option must have only 'c', 't', 'n' and/or '6'";
    }
    $do_weighting{c} = 1 if ($opt_weigh =~ /c/);
    $do_weighting{t} = 1 if ($opt_weigh =~ /t/);
    $do_weighting{n} = 1 if ($opt_weigh =~ /n/);
    $do_weighting{score} = 1 if ($opt_weigh =~ /6/);
}

#######################################################################

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}

#######################################################################

say STDERR "reading $opt_gff into GFF::Tree...";
my $gt = GFF::Tree->new(file => $opt_gff,normalize => 1, lenient => 1);
say STDERR "done";

#######################################################################
# create scores hash

my %methylation; # { locus_id => [seq, start, end, strand, c, t, n, score, lines] }

my $pp = GFF::Parser->new(file => $opt_gff);
while (defined(my $gff = $pp->next())){
    if (defined(my $locus = $gff->get_column($opt_tag))){
        $methylation{$locus}=[$gff->sequence,$gff->start,$gff->end,$gff->strand,0,0,0,0,0]; # c, t, n, score, lines
    }
}

#######################################################################
# Main loop

my $input_parser = GFF::Parser->new(file => $opt_input,normalize => 1);
my $counter = Counter->new();

while (defined(my $gff = $input_parser->next())){
    my ($seq, $start, $end, $strand, $score, $c, $t, $n) 
    = ($gff->sequence, $gff->start, $gff->end, $gff->strand, $gff->score(), $gff->get_column('c'), $gff->get_column('t'),$gff->get_column('n'));
    $strand //= '.';
    my $length = $end - $start + 1;

    my $gffstring = $gff->to_string;

    # if ($c+$t==0){ die "$gffstring"; }
    my @results = $gt->search_overlap($seq,$start,$end);

    for my $result (@results){
        my $overlap = $result->{overlap};
        my $weighed_c = defined $c && $do_weighting{c} ? $c * $overlap / $length : $c;
        my $weighed_t = defined $t && $do_weighting{t} ? $t * $overlap / $length : $t;
        my $weighed_n = defined $n && $do_weighting{n} ? $n * $overlap / $length : $n;
        my $weighed_score = defined $score && $do_weighting{score} ? $score * $overlap / $length : $score;

        if (defined(my $overlapped_locus = $result->{item}->get_column($opt_tag))){
            # locus may not have been in present in original annotation, so we'll add it here.
            if (! exists $methylation{$overlapped_locus}){
                $methylation{$overlapped_locus}=[$seq,$start,$end,$strand,0,0,0,0,0];
            }
            $methylation{$overlapped_locus}[4]+=$weighed_c if (defined $weighed_c);
            $methylation{$overlapped_locus}[5]+=$weighed_t if (defined $weighed_t);
            $methylation{$overlapped_locus}[6]+=$weighed_n if (defined $weighed_n);
            $methylation{$overlapped_locus}[7]+=$weighed_score if (defined $weighed_score);
            $methylation{$overlapped_locus}[8]++;
        }
        else{
            die sprintf("%s does not have an $opt_tag field?", $result->{item}->to_string);
        }
    }
    $counter->increment();
}

#######################################################################
# output

sub fmt{ $_[0] =~ /^\d+$/ ? $_[0] : sprintf("%.4f", $_[0]); }

LOCUS:
for my $id (sort keys %methylation) {
    my ($seq, $start, $end, $strand, $c, $t, $n, $score, $lines) = @{$methylation{$id}};
    my $length = $end - $start + 1;

    next LOCUS if ($lines == 0 && ! $opt_no_skip);

    my @attributes = ("ID=$id");
    if ($c > 0 || $t > 0){
        push @attributes, "c=" . fmt($c), "t=" . fmt($t);
    }
    if ($n > 0){
        push @attributes, "n=" . fmt($n);
    }
    # print score if non-zero or reporting derivative score
    if (($score > 0 || $opt_report_in_score =~ /^score_/) && $opt_report_in_score ne '6' && $opt_report_in_score ne 'score'){
        push @attributes, "score=" . fmt($score);
    }
    my $attr_string = join ";", @attributes;

    my $report_score;
    given ($opt_report_in_score){
        when ('c') { $report_score = $c; }
        when ('t') { $report_score = $t; }
        when ('n') { $report_score = $n; }
        when ('methyl') { $report_score = $c + $t > 0 ? $c / ($c + $t) : 0; }
        when ($_ eq '6' || $_ eq 'score') { $report_score = $score; }
        when ($_ eq 'l' || $_ eq 'lines') { $report_score = $lines; }
        when ('score_over_n')      { $report_score = $n > 0 ? $score / $n : 0; }
        when ('score_over_lines')  { $report_score = $lines > 0 ? $score / $lines : 0; }
        when ('score_over_length') { $report_score = $length > 0 ? $score / $length : 0; }
        default { die "unknown --report-in-score value $_" }
    }

    $strand //= '.';
    say join "\t", $seq, 'win', $opt_feature, $start, $end, fmt($report_score), $strand, $lines, $attr_string;
}

=head1 NAME

window_by_annotation_2.pl - a flexible window-by-annotaion.

=head1 What it does:

 For each input GFF entry G:
   For each annotation GFF entry A overlapping G:
     Record 'c', 't', 'n', and the 'score' values of G in A.
     Weigh by overlap percentage as specified by --weigh option (see details below)
  
 For each annotation GFF entry A:
   Report a GFF line containing:
     The locus tag, 'c', 't', 'n' in column 9.
     The number of lines in column 8
     The --report-in-score value in column 6 (see details below)

=head1 Usage examples:

Ordinary methylation windowing on single-c input:

 window_by_annotation_2.pl -g anno.gff -k -t 'ID' -o windowed.gff single-c.gff

Methylation windowing on windowed input.  Note that c, t are weighed:

 window_by_annotation_2.pl -g anno.gff -k -t 'ID' -w ct -o windowed.gff w50.gff

RNA-seq counting. Note that we are reporting the line-count in scores with '-s l':

 window_by_annotation_2.pl -g anno.gff -k -t 'ID' -s l -o windowed.gff w50.gff

Average score per annotation window, average by total 'n':

 window_by_annotation_2.pl -g anno.gff -k -t 'ID' -s score_over_n -o windowed.gff single-c.gff

Average score per annotation window, average by line count:

 window_by_annotation_2.pl -g anno.gff -k -t 'ID' -s score_over_lines -o windowed.gff single-c.gff

Average score per annotation window, average by annotation window length:

 window_by_annotation_2.pl -g anno.gff -k -t 'ID' -s score_over_length -o windowed.gff single-c.gff

Average weighed score per annotation window (useful for windowed input), average by line count:

 window_by_annotation_2.pl -g anno.gff -k -t 'ID' -s score_over_lines -w s -o windowed.gff single-c.gff

=head1 OPTIONS

=over 2 


=item  -g <file> | --gff <file>

Annotation file.

=for Euclid
    file.type:        readable

=item  <input>

=for Euclid
    input.type:        readable

=item  -t <tag> | --tag <tag>

Locus tag in --gff annotation file. Defaults to 'ID'.

=item -o <file> | --output <file>

=for Euclid
    file.default:     '-'

=item  -f <featurename> | --feature <featurename>

=for Euclid
    featurename.default:     'window'

=item  -s <value> | --report-in-score <value>

What should be reported in scores column? Default 'methyl'. Valid values are: 

=for Euclid
    value.default:     'methyl'

=over 4

=item c (The sum of 'c')

=item t (The sum of 't')

=item n (The sum of 'n')

=item l (Number of overlapping GFF entries/lines)

=item lines (Same as above)

=item 6 (The sum of scores (column 6))

=item score (Same as above)

=item methyl (c / (c + t))

=item score_over_n (score divided by 'n')

=item score_over_lines (score divided by number of lines)

=item score_over_length (score divided by length of annotation)

=back

=item  -k | --no-skip 

=item  -w <values> | --weigh <values>

Weigh the specified values by the amount of overlap by annotation window.  '6'
means weigh the scores columns, 'n', 'c', and 't' means weigh n, c, and t,
respectively. Can be combined ('6ctn' means weigh everything, '6n' means weigh
score and 'n').

For example, if you pass '-w s' and an input window-50 GFF entry is overlap by
an annotation entry by 20 bases, 40% of the score of the input would contribute
to the annotation.

=item --help | -h

=back

=cut

