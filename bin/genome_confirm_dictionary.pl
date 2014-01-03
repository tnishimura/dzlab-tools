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
use lib "$FindBin::Bin/../lib";
use FastaReader;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || ! $opt_from || !$opt_to || !$opt_dictionary;

my $config = LoadFile($opt_dictionary);
my %alignment = %{$config->{alignment}};

my $from = FastaReader->new(file => $opt_from, slurp => 1);
my $to   = FastaReader->new(file => $opt_to,   slurp => 1);

for my $chr (sort keys %alignment) {
    my $from_length = $from->get_length($chr);
    my $to_length   = $to->get_length($chr);
    my $from_coverage = 0;
    my $to_coverage = 0;
    my $diff_total = 0;

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
        $from_coverage += $from_end - $from_start + 1;
        $to_coverage   += $to_end - $to_start + 1;
        $diff_total += $diff;

        # if ($diff > 100){
        #     die 
        #     $from->get($chr, $from_start, $from_end) .
        #     " ========= \n" .
        #     $to->get($chr, $to_start, $to_end)
        # }
    }
    printf("summary: %s from coverage: %d / %d = %.2f\n", $chr, $from_coverage, $from_length, $from_coverage / $from_length);
    printf("summary: %s to coverage: %d / %d = %.2f\n", $chr, $to_coverage, $to_length, $to_coverage / $to_length);
    say "summary: total # of diffs: $diff_total";
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

Conversion dictionary, output by genome_make_dictionary.pl These are
*.alignment files found in the coord/ directory in the install path
(c:\dzlab-tools\coord on windows).  

=item -h | --help

=back

=cut

