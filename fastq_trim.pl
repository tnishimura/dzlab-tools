#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use English;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if ! ($opt_input || ! $opt_adapter );

$opt_header_length //= length($opt_adapter);

my $header_seq = substr $opt_adapter, 0, $opt_header_length;
my $header = q/$header_seq/;

my ($counter, $too_short, $too_long, $no_header) = (0,0,0,0);

open my $fh, '<', $opt_input;
my $output_file = $opt_out || $opt_input . '.trimmed';
open my $outfh, '>', $output_file;


while (defined (my $line = readfastq($fh))){
    INFO("$counter") if ++$counter % 100000 == 0;
    if ($line->[1] =~ m/$header_seq/){
        my $p = $PREMATCH;
        my $len = length $p;

        if (defined $opt_minimum_insert_size && $len < $opt_minimum_insert_size){
            ++$too_short;
        }
        elsif (defined $opt_maximum_insert_size && $len > $opt_maximum_insert_size){
            ++$too_long;
        } else{
            printf $outfh "%s\n%s\n%s\n%s\n", $line->[0], $p, $line->[2],  (substr $line->[3], 0, $len);
        }
    } else{
        ++$no_header;
    }
}
close $fh;
close $outfh;

INFO("total number of read: $counter");
INFO("inserts that were too short: $too_short");
INFO("inserts that were too long:  $too_long");
INFO("inserts that did not have the adapter header: $no_header");
INFO("number of passing reads: " . ($counter - $too_short - $too_long - $no_header));
INFO("percentage of passing reads: " . ($counter - $too_short - $too_long - $no_header) / $counter);

sub readfastq{
    my $fh = shift;
    my @line = map { my $x = scalar <$fh>; $x =~ tr/\r\n//d if defined $x; $x } (1..4);
    if (grep { defined } @line){
        return \@line;
    }
    return;
}

sub INFO{
    my $msg = shift;
    say STDERR $msg;
}

=head1 NAME

fastq_trim.pl - Trim out adapter for fastq reads.

=head1 SYNOPSIS

trim off adapter by finding 10 bp header from default adapter sequence and cutting it and trailing bases off:

 fastq_trim.pl -i input.fastq -h 10 

same as above but use a custom adapter sequence:

 fastq_trim.pl -i input.fastq -h 10 -a ACCCGGTGAGACATGAC

Even when you pass a custom adapter, only '-h' bases will be used.

=head1 OPTIONS

=over

=item  --adapter <seq> | -a <seq>

Adapter sequence.  For example:

 -a AGATCGGAAGA

=item --header-length <len> | -h <len>

Number of bases from the adapter sequence to use to search.  Default to entire adapter length.

=for Euclid
    len.type: int > 0 

=item  --minimum-insert-size <size> | -m <size>

Don't report post-trim reads under this long. 

=for Euclid
    size.type:    int > 0

=item  --maximum-insert-size <size> | -M <size>

Don't report post-trim reads over this long. 

=for Euclid
    size.type:    int > 0

=item --input <file> | -i <file>

Input file

=item --out <file> | -o <file>

Output file.  Default to <input_file>.trimmed

=back

=cut

