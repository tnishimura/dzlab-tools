#!/usr/bin/env perl
use 5.10.0;
use strict;
use warnings FATAL => "all";
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;

END {close STDOUT}
$|=1;


say "file\tcount\tmean\tmedian";

sub sum{
    my $arr = shift;
    my $total = 0;
    my $count = 0;
    for (@$arr) {
        if (defined $_){
            $total += $_;
            $count += 1;
        }
    }
    return ($total, $count);
}


for my $file (@ARGV) {
    my $filesize = (stat $file)[7];
    my $index = 0;
    my $p = GFF::Parser->new(file => $file);
    open my $fh, '<', $file;
    my @cplust;
    $#cplust = int $filesize/50;
    #while (defined(my $gff = $p->next())){
    while (defined(my $line = <$fh>)){
        chomp $line;
        if ($line =~ /c=(\d+);t=(\d+)/){
            $cplust[$index++] = $1+$2;
        }
    }
    @cplust = sort { $a <=> $b } @cplust;
    #say Dumper \@cplust;
    my ($sum, $length) = sum(\@cplust);
    my $mean   = $length ? $sum / $length : 0;
    my $median = $length ? $cplust[int($length/2)]: 0;
    printf "%s\t%d\t%.4f\t%.4f\n", $file, $length, $mean, $median;
    @cplust = ();
}

