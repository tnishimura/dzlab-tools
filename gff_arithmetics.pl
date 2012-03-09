#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use Getopt::Long;
use Scalar::Util qw/looks_like_number/;
use Pod::Usage;

use FindBin;
use lib "$FindBin::Bin/lib";
use Launch qw/cast/;

use File::Temp qw/tempfile/;

END {close STDOUT}
$| = 1;

my $result = GetOptions (
    "propagate-dots|d" => \(my $propagate_dots),
    "add-attributes|a" => \(my $add_attributes),
    "already-sorted|s" => \(my $assume_sorted),
);

my $expression = shift;
my $file_a = shift;
my $file_b = shift;

pod2usage(-verbose => 99) unless ($result && defined $expression && defined $file_a && defined $file_b);

if (! $assume_sorted){
    for ($file_a, $file_b){
        my $tmp = [tempfile("$_-XXXXXX", UNLINK => 1)]->[1];
        cast("sort -k1,1 -k4,4n -i $_ -o $tmp");
        rename $tmp, $_;
    }
}

#######################################################################
# main body

my $parser_left = GFF::Parser->new(file => $file_a, normalize => 0);
my $parser_right = GFF::Parser->new(file => $file_b, normalize => 0);

LEFT:
while (defined(my $left = $parser_left->next)){
    RIGHT:
    while (defined(my $right = $parser_right->next)){
        if(GFF::start_position_equal($left,$right)){

            say join "\t", 
            $left->sequence(), 
            $left->source() // '.', 
            $left->feature() // '.', 
            $left->start(), 
            $left->end(), 
            eval_expression($left->score(), $right->score()), 
            $left->strand() // '.', 
            $left->frame() //'.',
            combine_attributes($left, $right);
            #"asdf";
            last RIGHT;
        }
        elsif(GFF::start_position_lessthan($left,$right)){
            $parser_right->rewind($right);
            say $left->to_string;
            last RIGHT;
        }
        elsif(GFF::start_position_lessthan($right,$left)){
            say $right->to_string;
        }
    }
}
while (defined(my $right = $parser_right->next)){
    $right->score(-$right->score);
    say $right->to_string;
}

#######################################################################
# score combiner

sub eval_expression{
    my ($left_score, $right_score) = @_;
    no warnings;

    if ((! defined $left_score || !defined $right_score ) && $propagate_dots){
        return '.';
    }
    
    local $a = $left_score // 0;
    local $b = $right_score // 0;
    my $result = eval $expression;
    if ($@){
        if ($@ =~ /Illegal division by zero/){
            return '.';
        }
        else{
            die "expression error: $@";
        }
    }
    else{
        return $result;
    }
}

#######################################################################
# attribute combiner


sub combine_attributes{
    my ($left_gff, $right_gff) = @_;

    my @left_keys = $left_gff->list_attribute();
    my @right_keys = $right_gff->list_attribute();

    # find keys in both
    my %common_keys = do {
        my %key_counter;
        for (@left_keys, @right_keys){
            $key_counter{$_}++;
        }
        map{ $_ => 1 } grep { $key_counter{$_} > 1 } keys %key_counter;
    };

    # collect keys with identical values in both
    my %identical_keys = do{
        map {
            $_ => $left_gff->get_attribute($_) 
        }
        grep {
            ! (
                looks_like_number($left_gff->get_attribute($_)) && 
                looks_like_number($right_gff->get_attribute($_))
            ) 
            && $left_gff->get_attribute($_) eq $right_gff->get_attribute($_)
        }
        keys %common_keys;
    };

    my @accum;

    # print identical only once each
    for my $k (sort keys %identical_keys) {
        push @accum, "$k=" . $identical_keys{$k};
    }

    if ($add_attributes){
        # collect keys with numeric values in both, and add them
        my %added_keys = do{
            map {
                $_ => ($left_gff->get_attribute($_) + $right_gff->get_attribute($_))
            }
            grep {
                looks_like_number($left_gff->get_attribute($_)) && 
                looks_like_number($right_gff->get_attribute($_))
            }
            keys %common_keys;
        };

        # print left-only keys
        for my $k (@left_keys) {
            if (! exists $identical_keys{$k} && ! exists $added_keys{$k}){
                push @accum, "$k=" . $left_gff->get_attribute($k);
            }
        }
        # print right-only keys
        for my $k (@right_keys) {
            if (! exists $identical_keys{$k} && ! exists $added_keys{$k}){
                push @accum, "$k=" . $right_gff->get_attribute($k);
            }
        }
        # print added keys
        for my $k (sort keys %added_keys) {
            push @accum, "$k=" . $added_keys{$k};
        }
        return join ";", @accum;
    }
    else{
        # print left-only keys
        for my $k (@left_keys) {
            if (! exists $identical_keys{$k}){
                push @accum, "$k=" . $left_gff->get_attribute($k);
            }
        }
        # print right-only keys
        for my $k (@right_keys) {
            if (! exists $identical_keys{$k}){
                push @accum, "$k=" . $right_gff->get_attribute($k);
            }
        }
        return join ";", @accum;
    }
}


=head1 NAME

gff_arithmetics.pl - perform base-by-base arithmetic on scores column

=head1 SYNOPSIS

Usage: 
 
 gff_arithmetics.pl "$a + $b" file1.gff file2.gff

This script looks at two gff files line by line, and for each common base
(meaning the sequence and start position are the same), performs the arithmetic
expression, and prints a combined line. Any lines which are in one but not the
other are output as-is.

For example, say file1.gff contains:

 Chr1	TAIR8	gene	1	50	1	-	.	ID=first;c=1
 Chr1	TAIR8	gene	101	150	3	+	.	ID=third;c=3
 Chr1	TAIR8	gene	151	200	4	+	.	ID=fourth;c=123
 Chr1	TAIR8	gene	201	250	5	-	.	ID=fifth
 Chr1	TAIR8	gene	251	300	6	-	.	ID=sixth
 Chr1	TAIR8	gene	301	350	.	-	.	ID=seventh
 Chr1	TAIR8	gene	351	400	.	-	.	ID=eigth

And file2.gff contains:

 Chr1	TAIR8	gene	1	50	0	-	.	ID=first;c=92
 Chr1	TAIR8	gene	51	100	1	-	.	ID=second;c=7
 Chr1	TAIR8	gene	101	150	2	+	.	ID=third;c=45
 Chr1	TAIR8	gene	151	200	2	+	.	ID=fourth;c=3
 Chr1	TAIR8	gene	251	300	3	-	.	ID=sixth
 Chr1	TAIR8	gene	301	350	.	-	.	ID=seventh
 Chr1	TAIR8	gene	351	400	.	-	.	ID=eigth

Notice that there are positions unique to both.  Then,

 gff_arithmetics.pl '$a + $b' file1.gff file2.gff 

produces:

 Chr1    TAIR8   gene    1       50      1       -       .       ID=first;c=1;c=92
 Chr1    TAIR8   gene    51      100     1       -       .       ID=second;c=7
 Chr1    TAIR8   gene    101     150     5       +       .       ID=third;c=3;c=45
 Chr1    TAIR8   gene    151     200     6       +       .       ID=fourth;c=123;c=3
 Chr1    TAIR8   gene    201     250     5       -       .       ID=fifth
 Chr1    TAIR8   gene    251     300     9       -       .       ID=sixth
 Chr1    TAIR8   gene    301     350     0       -       .       ID=seventh
 Chr1    TAIR8   gene    351     400     0       -       .       ID=eigth

 gff_arithmetics.pl -a '$a * $b' file1.gff file2.gff 

produces:

 Chr1    TAIR8   gene    1       50      0       -       .       ID=first;c=93
 Chr1    TAIR8   gene    51      100     1       -       .       ID=second;c=7
 Chr1    TAIR8   gene    101     150     6       +       .       ID=third;c=48
 Chr1    TAIR8   gene    151     200     8       +       .       ID=fourth;c=126
 Chr1    TAIR8   gene    201     250     5       -       .       ID=fifth
 Chr1    TAIR8   gene    251     300     18      -       .       ID=sixth
 Chr1    TAIR8   gene    301     350     0       -       .       ID=seventh
 Chr1    TAIR8   gene    351     400     0       -       .       ID=eigth

Notice that the -a option adds the 'c's in column 9.

=head1 EXPRESSION

The expression can be any legal perl statement with $a holding the score of the
first file, and $b holding the score of the second.  In general, most
arithmetic statements that you might type into a graphing calculator will work:

 $a + $b               # addition
 $a * $b               # multiplication
 ($a * $b) * 3 + $a/$b # compound statement

Note: 'dots' are interpretted as zero, unless you pass the -d (or
--propagate-dots) option, in which case any position with a dot in one file is
output as another dot.

Also note: division by zero are detected and replaced with a dot.

=head1 OPTIONS

=over

=item -d | --propagate-dots

By default, dots (the GFF way of saying 'undefined') are interpretted as zero.
This may be desirable if you are, say, adding or subtracting.  But with this
option, any position with dots in either side are output as dots.

=item -a | --add-attributes

If there are common numeric attributes in both files, add them.  Usefull when, say,

=item  -s | --already-sorted 

By default, the script will sort both files. If you know that they are, you can
pass this option to skip this step.  Warning: if they aren't sorted, this
script will produce garbage, silently.

=back

=cut

