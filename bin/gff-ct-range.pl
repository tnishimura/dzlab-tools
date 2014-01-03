#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
END {close STDOUT}
$| = 1;
use FindBin;
use lib "$FindBin::Bin/../lib";
use GFF::Parser;
use Getopt::Long;

my @nuc;
my @mit;
my @chr;

my $result = GetOptions (
    "nuc|n=s{2}" => \@nuc,
    "mit|m=s{2}" => \@mit,
    "chr|c=s{2}" => \@chr,
);

if (!$result or !@nuc or !@mit or !@chr) {
    say "$0 --nuc <min> <max> --mit <min> <max> --chr <min> <max> <file>";
    exit 1;
}

my $parser = GFF::Parser->new(file => \*ARGV);

sub between{
    my ($start, $val, $end) = @_;
    return (
        ($start == -1 || $start <= $val) && 
        ($end   == -1 || $val <= $end)
    );
}

my $counter = 0;
while (defined(my $gff = $parser->next())){
    say STDERR $counter if (++$counter % 50000 == 0);
    my $seq = $gff->sequence();
    my $c = $gff->get_attribute('c');
    my $t = $gff->get_attribute('t');
    my $print = 0;
    if (defined $c and defined $t){
        my $ct = $c + $t;
        if ($seq =~ /chrc/i && between($chr[0], $ct, $chr[1])){
            say $gff;
        }
        elsif ($seq =~ /chrm/i && between($mit[0], $ct, $mit[1])){
            say $gff;
        }
        elsif (between($nuc[0], $ct, $nuc[1])){
            say $gff;
        }
    }
}

=head1 NAME

gff-iqr.pl - Read in GFF, print out lines within a certain C+T range

=head1 SYNOPSIS

Keep only gff lines where nuclear c+t is 10 to 20, mitochondrial is 100 to 200,
and chloroplast is 1 to 20.

 gff-iqr.pl --nuc 10 20 --mit 100 200 --chr 1 20 in.gff > out.gff

Keep only gff lines where nuclear c+t is at least 10, mitochondrial is at most 200,
and any chloroplast.

 gff-iqr.pl --nuc 10 -1 --mit -1 200 --chr -1 -1 in.gff > out.gff

=cut

