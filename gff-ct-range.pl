#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
END {close STDOUT}
$| = 1;
use FindBin;
use lib "$FindBin::Bin/lib";
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

gff-iqr.pl - 

=head1 SYNOPSIS

Usage examples:

 gff-iqr.pl [options]...

=head1 REQUIRED ARGUMENTS

=over

=item  --iqr <p25> <p75>

=back

=cut

