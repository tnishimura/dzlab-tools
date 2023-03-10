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
use List::Util qw/max min/;
use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "output|o=s"         => \(my $output),
    "input|i=s"          => \(my $input),
    "window-size|w=i"    => \(my $window_size),
    "distribute-score|d" => \(my $distribute_score),
    "always-align|a"     => \(my $always_align),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || (@ARGV == 0 && ! defined $input && -t STDIN)) ;  
my $p = GFF::Parser->new(normalize => 0, file => defined($input) ? $input : \*ARGV);

if (defined $output){
    open my $outfh, '>', $output;
    select $outfh;
}

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
        $gff->attribute_string // '.'
    ); 
    my $start = $gff->start();
    my $end = $gff->end();
    my $length = $end - $start + 1;
    if (looks_like_number $gff->start && looks_like_number $gff->end && $gff->start < $gff->end){
        if ($window_size){
            my $aligned_start = $start - ( $start % $window_size ) + 1; 
            my $aligned_end   = $aligned_start + $window_size - 1;
            while ($aligned_start <= $end){
                my $window_start = $always_align ? $aligned_start : max($start, $aligned_start);
                my $window_end   = $always_align ? $aligned_end : min($end, $aligned_end);

                my $score = $gff_split[5];
                if ($distribute_score){
                    $score = $score / $length * (min($end, $aligned_end) - max($start, $aligned_start) + 1);
                    if ($score !~ /^\d+$/){
                        $score = sprintf("%.6f", $score);
                    }
                }
                @gff_split[3,4] = ($window_start, $window_end);
                say join "\t", @gff_split[0 .. 4], $score, @gff_split[6 .. 8];

                $aligned_start += $window_size;
                $aligned_end += $window_size;
            }
        }
        else{
            if ($distribute_score){
                $gff_split[5] /= ($gff->end - $gff->start + 1);
            }
            for my $index ($gff->start .. $gff->end){
                @gff_split[3,4] = ($index, $index);
                say join "\t", @gff_split;
            }
        }
    }
}

if (defined $output){
    close \*STDOUT;
}

=head1 NAME

gff_explode.pl - explode windows into single-nucleotide or multi-nucleotide windows

=head1 SYNOPSIS

Explode into 1 bp:

 gff_explode.pl input.gff > output.gff 
 gff_explode.pl -o output.gff input.gff 
 gff_explode.pl -o output.gff -i input.gff 

For example the above would explode, 
 
 Chr1    TAIR8   gene    3631    5899    .       +       .       ID=AT1G01010

into: 

 Chr1    TAIR8   gene    3631    3631    1       +       .       ID=AT1G01010
 Chr1    TAIR8   gene    3632    3632    1       +       .       ID=AT1G01010
 Chr1    TAIR8   gene    3633    3633    1       +       .       ID=AT1G01010
 <lots of lines omitted>
 Chr1    TAIR8   gene    5897    5897    1       +       .       ID=AT1G01010
 Chr1    TAIR8   gene    5898    5898    1       +       .       ID=AT1G01010
 Chr1    TAIR8   gene    5899    5899    1       +       .       ID=AT1G01010

Explode into 50 bp:

 gff_explode.pl -w 50 input.gff        

For example the above would explode, 
 
 Chr1    TAIR8   gene    3631    5899    .       +       .       ID=AT1G01010

into: 

 Chr1    TAIR8   gene    3631    3650    1       +       .       ID=AT1G01010
 Chr1    TAIR8   gene    3651    3700    1       +       .       ID=AT1G01010
 Chr1    TAIR8   gene    3701    3750    1       +       .       ID=AT1G01010
 <lots of lines omitted>
 Chr1    TAIR8   gene    5751    5800    1       +       .       ID=AT1G01010
 Chr1    TAIR8   gene    5801    5850    1       +       .       ID=AT1G01010
 Chr1    TAIR8   gene    5851    5899    1       +       .       ID=AT1G01010

Explode into 1 bp, distributing scores amongst exploded fragments

 gff_explode.pl -d input.gff     

Explode into 50 bp, distributing scores amongst exploded fragments

 gff_explode.pl -w 50 -d input.gff     


=head1 OPTIONS

=over

=item  -o <file> | --output <file>

=item -w <width> | --window-size <width>

Explode into multi-bp windows instead of single-bp positions.  windows are
always aligned at even multiples of window-size.  For example, an input GFF
entry from 1234 to 1333 with window size 50 will result in windows for
1234-1250, 1251-1300, 1301-1350.

=item  -d | --distribute-score

If enabled, score will be distributed amongst the exploded fragments (for
example, if a length 123 window with score 10 is exploded into 123 single-bp
positions, normally each position would retain a score of 10.  With this
option, each position would only get 10/123).

=item -a | --always-align

If enabled, always align windows to the window_size. For example, 

gff_explode.pl -w 50 -a input.gff        

 Chr1    TAIR8   gene    3601    3650    1       +       .       ID=AT1G01010
 Chr1    TAIR8   gene    3651    3700    1       +       .       ID=AT1G01010
 Chr1    TAIR8   gene    3701    3750    1       +       .       ID=AT1G01010
 <omitted...>
 Chr1    TAIR8   gene    5751    5800    1       +       .       ID=AT1G01010
 Chr1    TAIR8   gene    5801    5850    1       +       .       ID=AT1G01010
 Chr1    TAIR8   gene    5851    5900    1       +       .       ID=AT1G01010

Note that score distribution, if enabled, always distributes according to the
real start and end of the input gff, not the aligned start/end.

=item --help | -h

=back

=cut
