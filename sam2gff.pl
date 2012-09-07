#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;

use FindBin;
use lib "$FindBin::Bin/lib";
use SamParser;

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "print-unmapped|u" => \(my $print_unmapped),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || (! @ARGV && -t STDIN));  

my %stats;

while (defined(my $line = <ARGV>)){
    chomp $line;
    my $r = parse_sam_line($line);

    next if (ref $r ne 'SAMLINE');

    my $mapped     =  $r->{mapped}     ;
    my $seqid      =  $r->{seqid}      ;
    my $length     =  $r->{length}     ;
    my $leftmost   =  $r->{leftmost}   ;
    my $rightmost  =  $r->{rightmost}  ;
    my $failed_qc  =  $r->{failed_qc}  ;
    my $reverse    =  $r->{reverse}    ;
    my $readid     =  $r->{readid}     ;
    my $mapq       =  $r->{mapq}     ;
    my $seq        =  $r->{seq}        ;
    my $cigar      =  $r->{cigar}      ;

    next if (! $mapped && ! $print_unmapped);

    if ($mapped){
        $stats{$seqid}++;
        say(join "\t",
            $seqid,
            $length,
            $readid,
            $leftmost,
            $rightmost, 
            $mapq,
            ($reverse ? '-' : '+'),
            q{.},
            "cigar=$cigar;seq=$seq"
        );
    }
    else{
        say(join "\t",
            $seqid,
            q{.},
            q{.},
            q{.},
            q{.},
            q{.},
            q{.},
            q{.},
            "seq=$seq"
        );
    }
}

for my $seq (sort keys %stats) {
    say STDERR "$seq: $stats{$seq}";
}

=head1 NAME

sam2gff.pl in.sam > out.gff

=cut

