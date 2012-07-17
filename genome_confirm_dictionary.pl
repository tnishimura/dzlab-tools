#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use YAML qw/LoadFile DumpFile/;
use FindBin;
use lib "$FindBin::Bin/lib";
use FastaReader;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || ! $opt_from || !$opt_to || !$opt_dictionary;

my $config = LoadFile("$FindBin::Bin/conf/$opt_dictionary.alignment");
my %alignment = %{$config->{alignment}};

my $from = FastaReader->new(file => $opt_from, slurp => 1);
my $to   = FastaReader->new(file => $opt_to,   slurp => 1);

for my $chr (sort keys %alignment) {
    for my $fragment (@{$alignment{$chr}}){
        my ($from_start, $from_end, $to_start, $to_end) = @{$fragment};
        my $size = abs($from_end - $from_start + 1);
        die "unequal size $chr($from_start, $from_end, $to_start, $to_end)" unless $from_start - $from_end == $to_start - $to_end;

        my $diff = 
        change(
            $from->get($chr, $from_start, $from_end), 
            $to->get($chr, $to_start, $to_end)
        );

        say "$chr ($from_start, $from_end, $to_start, $to_end) = $diff / $size";

        # if ($diff > 100){
        #     die 
        #     $from->get($chr, $from_start, $from_end) .
        #     " ========= \n" .
        #     $to->get($chr, $to_start, $to_end)
        # }
    }
}

sub change{
    my ($x, $y) = @_;
    my @x_split = split //, $x;
    my @y_split = split //, $y;

    my $diff = 0;
    for my $i (0 .. $#x_split) {
        if ($x_split[$i] ne $y_split[$i]){
            $diff++;
        }
    }
    return $diff;
}


=head1 NAME
 
gff_convert_coord.pl - convert coordinates according to a dictionary.
 
=head1 SYNOPSIS

Convert column 4 and 5 in input.gff to output.gff according to rice-5.0-to-6.1:

 gff_convert_coord.pl -d rice-5.0-to-6.1 from.fasta to.fasta

=head1 OPTIONS

=over

=item --from <from>

=item --to <to>

=item  -d <dict> | --dictionary <dict>

Conversion dictionary.  These are files found in the coord/ directory in the
install path (c:\dzlab-tools\coord on windows).  The '.alignment' prefix is
optional. Currently available are:

 rice-5.0-to-6.1
 rice-6.1-to-5.0
 at-tair8-to-tair10
 at-tair10-to-tair8

=item -h | --help

=back

=cut

