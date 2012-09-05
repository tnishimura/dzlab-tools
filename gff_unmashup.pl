#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "input|i=s" => \(my $input_file),
    "window-size|w=i" => \(my $window_size = 1),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || ! $input_file);  

my $file_template = $input_file =~ s/(\.\w+)$//r . ".%s" . ($1 // '');

# read
open my $input, '<:crlf', $input_file;


# read header line
my @headers = do{
    my $header_line = <$input>;
    chomp $header_line;
    split /\t/, $header_line;
};

# create file handle for each
my %file_handles = map {
    my $header = $headers[$_];
    my $filename = sprintf($file_template, $header);
    open my $fh, '>:crlf', $filename;
    ($_ => $fh);
} (2 .. $#headers);

while (defined(my $line = <$input>)){
    chomp $line;
    my @parts = split /\t/, $line;
    my ($seq, $coord) = @parts;

    die "line $.: # of columns does not match header" if @parts != @headers;

    for (2 .. $#headers) {
        say {$file_handles{$_}} join "\t", $seq, qw/. ./, $coord, $coord + $window_size - 1, $parts[$_];
    }
}
close $input;

close $_ for values %file_handles;

=head1 gff_unmashup.pl 

Given an input mashup file, undo the mashup process by creating a gff file for
each column with the value as column 6.

Usage examples:

 gff_unmashup.pl -w 50 -i w50mashup.txt
 gff_unmashup.pl -w 1 -i w1mashup.txt

=over

=item  -i <file> | --input <file>

Input file.  No, you can pass it as STDIN b/c we need to know how to name the output files.

=item  -w <bp> | --window-size <bp>

Assume the windows are <bp> wide.  Default 1. 

=back

=cut

