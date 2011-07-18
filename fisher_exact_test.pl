#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
#use Math::BigFloat lib => 'GMP';
use List::Util qw/sum/;
use Pod::Usage;
use Getopt::Long;

my $sep = 'tab';
my $cola;
my $colb;
my $colc;
my $cold;
my $input;
my $output = q{-};
my $help;
my $pval = .000_1;
my $log_threshold = -200;
my $result = GetOptions (
    "input|i=s"     => \$input,
    "output|o=s"    => \$output,
    "column-a|a=i"  => \$cola,
    "column-b|b=i"  => \$colb,
    "column-c|c=i"  => \$colc,
    "column-d|d=i"  => \$cold,
    "p-value|p=i"   => \$pval,
    "threshold|t=i" => \$log_threshold,
    "sep|s=s"       => \$sep,
    "help"          => \$help,
);
pod2usage(-verbose => 99) 
if (!($cola && $colb && $colc && $cold && $input) || $log_threshold > 0 || $help || !$result);

# logfactorial(n) = log10(n!)
sub logfactorial {
    my ($n) = @_;
    return 0 if $n==0;
    my $accum = 0;
    my $log10 = log(10);
    for (my $i=1 ; $i<=$n ; $i++){
        $accum += log($i)/$log10;
    }
    return $accum;
};

sub fet_using_log{
    my ($a,$b,$c,$d) = @_;
    my $n = $a+$b+$c+$d;
    my $logp = sum (
        logfactorial($a+$b) , logfactorial($c+$d) , logfactorial($a+$c) , logfactorial($b+$d) , 
        -logfactorial($a)   , -logfactorial($b)   , -logfactorial($c)   , -logfactorial($d)   , -logfactorial($n) , 
    );
    return ($logp < $log_threshold ? 0 : sprintf("%g",10 ** $logp));
}

if ($sep eq 'tab'){
    $sep = "\t";
}

my $outfh;
if ($output eq q{-}){
    $outfh = *STDOUT;
}
else {
    open $outfh, '>', $output;
}

open my $fh, '<', $input;

my $counter = 1;
while (defined (my $line = <$fh>)){
    $line =~ tr/\n\r//d;
    my @parts = split /$sep/, $line;
    my ($a, $b, $c, $d) = @parts[$cola-1, $colb-1, $colc-1, $cold-1];

    say $outfh join $sep, @parts, fet_using_log($a,$b,$c,$d);
}

close $fh;
if ($output ne q{-}){
    close $outfh;
}

=head1 NAME

fisher_exact_test.pl - perform fisher exact test on file

=head1 SYNOPSIS

fisher_exact_test.pl -a 2 -b 4 -c 6 -d 8 -i input.txt -o output.txt

=head1 DESCRIPTION

 fisher_exact_test.pl 
 Input: 
 takes as input a multi-column file and a list of column numbers to identify
 a,b,c,d to calculate the p-value according the the Fisher Exact Test, illustrated below.  For each
 line in the in the in

 +-----+-----+-----+
 |  a  |  b  | a+b |
 +-----+-----+-----+
 |  c  |  d  | c+d |
 +-----+-----+-----+
 | a+c | b+d |  n  |
 +-----+-----+-----+
 
     (a+b)! * (c+d)! * (a+c)! * (b+d)!
 p = ---------------------------------
          a! * b! * c! * d! * n!

=head1 OPTIONS

 --column-a      -a  Position of column a in file.
 --column-b      -b  Position of column b in file.
 --column-c      -c  Position of column c in file.
 --column-d      -d  Position of column d in file.
 --sep           -s  File column separator (default: tab) 
 --significance  -p  If P-value is less than this, mark it with a '+' sign
                     Otherwise, '-'. Default: .0001
 --input         -i  Input file
 --output        -o  Output file (default to screen)
 --threshold     -t  Report 0 for p-value if it is less than 10^threshold.
                     (default -200).

 --help        print this information

=cut

