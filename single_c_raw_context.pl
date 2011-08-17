#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use FastaReader;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use File::Basename;
use List::Util qw/sum/;
END {close STDOUT}
$| = 1;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || ! defined $opt_reference || ! @opt_inputs;

my $fasta = FastaReader->new(file => $opt_reference, slurp => 1);

for my $file (@opt_inputs) {
    my $parser = GFF::Parser->new(file => $file, normalize => 0);

    (my $filename_template = $file) =~ s/\.gff$//;
    my %filehandles;
    my %counts;

    while (defined(my $gff = $parser->next())){
        my $forward = $fasta->get_context_raw($gff->sequence, $gff->start, rc => 0);
        my $reverse = $fasta->get_context_raw($gff->sequence, $gff->start, rc => 1);
        my $raw = ($forward =~ /^C/) ? $forward : $reverse;

        next unless (length $raw == 3);
        if (! exists $filehandles{$raw}){
            open $filehandles{$raw}, '>', "$filename_template-$raw.gff";
        }
        $gff->feature($gff->feature() . "-$raw");
        say {$filehandles{$raw}} $gff->to_string();
        ++$counts{$raw};
    }
    close $_ for values %filehandles;
    # log
    {
        open my $log, '>', "$filename_template-rawcont_log.txt";
        my $total = sum values %counts;
        for my $raw (sort keys %counts) {
            printf $log "%s\t%d\t%.4f\n", $raw, $counts{$raw}, $counts{$raw}/$total;
        }

        close $log;
    }
}


=head1 NAME

single_c_raw_context.pl - separate single-c file by raw triplet context.

=head1 SYNOPSIS

Usage examples:

 single_c_raw_context.pl -r genome.fasta 

=head1 REQUIRED ARGUMENTS

=over

=back

=head1 OPTIONS

=over

=item  -r <reference> | --reference <reference>

=for Euclid
    reference.type:        readable

=item <inputs>...

=for Euclid
    inputs.type:        readable

=item --help | -h

=back

=cut

