#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;

END {close STDOUT}
$| = 1;
use Pod::Usage;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../lib";
use DZUtil qw/safediv/;


my $result = GetOptions (
    "eland1|1=s" => \my $eland1,
    "eland2|2=s" => \my $eland2,
    "correl|c=s" => \my $correl,
    "tissue|t=s" => \my $tissue,
    "batch|b=s" => \my $batch,
);
pod2usage(-verbose => 1) if (!$eland1 || !$tissue || !$batch);

$batch  = 'batch-' . $batch;

my %elands = (left => $eland1, right => $eland2);
my %ecount = (
    left => {ALL => 0, NM => 0, U0 => 0, U1 => 0, U2 => 0, RT => 0, UT => 0},
    right => {ALL => 0, NM => 0, U0 => 0, U1 => 0, U2 => 0, RT => 0, UT => 0}
);

while (my ($side,$file) = each %elands) {
    next if ! defined $file;
    open my $fh, '<', $file;

    while (defined(my $line = <$fh>)){
        chomp $line;
        my @efields1 = split /\t/, $line;

        ++$ecount{$side}{ALL};

        ++$ecount{$side}{NM} if $efields1[2] =~ m/NM/;

        ++$ecount{$side}{U0} if $efields1[2] =~ m/^1:[09]:[09]/;
        ++$ecount{$side}{U1} if $efields1[2] =~ m/^[09]:1:[09]/;
        ++$ecount{$side}{U2} if $efields1[2] =~ m/^[09]:[09]:1/;
        ++$ecount{$side}{RT} if $efields1[2] !~ m/NM/ and $efields1[3] =~ m/,/;
    }
    close $fh;
    $ecount{$side}{UT} = $ecount{$side}{U0} + $ecount{$side}{U1} + $ecount{$side}{U2};
}

my %ccounts = (ALL => 0, NM => 0, UT => 0, RT => 0, RR => 0, UU => 0);
use Scalar::Util qw/looks_like_number/;

if ($correl) {
    open my $CORR, '<', $correl or die "Can't open $correl";
  CORREL_LINE:
    while (<$CORR>) {
        my @cfields = split /\t/, $_;
        next if @cfields != 9 || ! looks_like_number $cfields[5] ;


        ++$ccounts{ALL};

        ++$ccounts{NM} if $cfields[5]  < 1;
        ++$ccounts{UT} if $cfields[5] == 1;
        ++$ccounts{RT} if $cfields[5]  > 1;

        $ccounts{RR} += 1/2 if $cfields[1] =~ m{U\/R|R\/U} and $cfields[5] == 1;
        ++$ccounts{RR} if $cfields[1] =~ m{R\/R} and $cfields[5] == 1;
        ++$ccounts{UU} if $cfields[1] =~ m{U\/U} and $cfields[5] == 1;
    }
}


print join ("\t",
                 'TISSUE/BATCH',
                 '/1_TOT',
                 '/1_U0',
                 '/1_UT',
                 '/1_U0/UT',
                 '/2_TOT',
                 '/2_U0',
                 '/2_UT',
                 '/2_U0/UT',
                 'REPEAT_RESOLV',
                 'TOTAL_MATCHES',
                 'RR/TM'
             ), "\n" ;

print join ("\t",
                 join (q{|}, $tissue, $batch),
                 $ecount{left}{ALL},
                 $ecount{left}{U0},
                 $ecount{left}{UT},
                 safediv($ecount{left}{U0},$ecount{left}{UT}),
                 $ecount{right}{ALL},
                 $ecount{right}{U0},
                 $ecount{right}{UT},
                 safediv($ecount{right}{U0},$ecount{right}{UT}),
                 $ccounts{RR},
                 $ccounts{UT},
                 safediv($ccounts{RR}, $ccounts{UT}),
             ), "\n";

print "================================================================================\n";
print "$eland1\n";
print "================================================================================\n";
print "No. total reads:\t$ecount{left}{ALL}\n";
print "No. non-matched reads:\t$ecount{left}{NM}\n";
print "No. 0-mm matched reads:\t$ecount{left}{U0}\n";
print "No. 1-mm matched reads:\t$ecount{left}{U1}\n";
print "No. 2-mm matched reads:\t$ecount{left}{U2}\n";
print "No. matched reads:\t$ecount{left}{UT}\n";
print "No. repeat reads:\t$ecount{left}{RT}\n";
print "Ratio of matched reads in total reads:\t", safediv($ecount{left}{UT} , $ecount{left}{ALL}), "\n";
print "Ratio of 0-mm reads in total reads:\t", safediv($ecount{left}{U0} , $ecount{left}{ALL}), "\n";
print "Ratio of 0-mm reads in matched reads:\t", safediv($ecount{left}{U0} , $ecount{left}{UT}), "\n";
print "Ratio of non-matched reads in total reads:\t", safediv($ecount{left}{NM} , $ecount{left}{ALL}), "\n";
print "Ratio of repeats reads in total reads:\t", safediv($ecount{left}{RT} , $ecount{left}{ALL}), "\n";
print "================================================================================\n";

if (defined $elands{right}){
    print "================================================================================\n";
    print "$eland2\n";
    print "================================================================================\n";
    print "No. total reads:\t$ecount{right}{ALL}\n";
    print "No. non-matched reads:\t$ecount{right}{NM}\n";
    print "No. 0-mm matched reads:\t$ecount{right}{U0}\n";
    print "No. 2-mm matched reads:\t$ecount{right}{U1}\n";
    print "No. 2-mm matched reads:\t$ecount{right}{U2}\n";
    print "No. matched reads:\t$ecount{right}{UT}\n";
    print "No. repeat reads:\t$ecount{right}{RT}\n";
    print "Ratio of matched reads in total reads:\t", safediv($ecount{right}{UT} , $ecount{right}{ALL}), "\n";
    print "Ratio of 0-mm reads in total reads:\t", safediv($ecount{right}{U0} , $ecount{right}{ALL}), "\n";
    print "Ratio of 0-mm reads in matched reads:\t", safediv($ecount{right}{U0} , $ecount{right}{UT}), "\n";
    print "Ratio of non-matched reads in total reads:\t", safediv($ecount{right}{NM} , $ecount{right}{ALL}), "\n";
    print "Ratio of repeats reads in total reads:\t", safediv($ecount{right}{RT} , $ecount{right}{ALL}), "\n";
    print "================================================================================\n";
}

if ($correl) {
    print "================================================================================\n";
    print "$correl\n";
    print "================================================================================\n";
    print "No. total reads:\t$ccounts{ALL}\n";
    print "No. non-matched reads:\t$ccounts{NM}\n";
    print "No. matched reads:\t$ccounts{UT}\n";
    print "No. preserved unique reads:\t$ccounts{UU}\n";
    print "No. repeat reads:\t$ccounts{RT}\n"; #
    print "No. resolved repeats:\t$ccounts{RR}\n";
    print "Ratio of preserved unique reads in total matched reads:\t", safediv($ccounts{UU} , $ccounts{UT}), "\n";
    print "Ratio of matched reads in total reads:\t", safediv($ccounts{UT} , $ccounts{ALL}), "\n";
    print "Ratio of non-matched reads in total reads:\t", safediv($ccounts{NM} , $ccounts{ALL}), "\n";
    print "Ratio of repeat reads in total reads:\t", safediv($ccounts{RT} , $ccounts{ALL}), "\n"; #
    print "Ratio of resolved repeats in matched reads:\t", safediv($ccounts{RR} , $ccounts{UT}), "\n";
    print "================================================================================\n";
}

