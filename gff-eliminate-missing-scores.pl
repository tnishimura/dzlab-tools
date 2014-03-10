#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;
use List::MoreUtils qw/all/;
use Pod::Usage;
use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::DataFrame;

my $result = GetOptions (
    "suffix|s=s" => \(my $suffix = ".no-missing-scores.gff"),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || !@ARGV);  

my @files = @ARGV;

my %output = map 
{ 
    my $input_file = $_; 
    my $output_file = ($input_file =~ s/\.gff$//r ) . $suffix; 
    open my $fh, '>', $output_file; 
    $input_file => $fh 
} @files;

my $frame = GFF::DataFrame->new_from_gffs(@files);

my $iter = $frame->make_iterator;
while (my @results = $iter->()){
    my ($sequence, $start, $end, %files2gffs) = @results;
    if (all { defined($files2gffs{$_}) and defined $files2gffs{$_}->score } @files){
        for my $f (@files) {
            say {$output{$f}} $files2gffs{$f};
        }
    }
}

close $_ for values %output;

=head1 gff-eliminate-missing-scores.pl 

Given a set of gff's, create another gff for each file which keeps only positions which have scores in all files.
The following call:

 gff-eliminate-missing-scores.pl file1.gff file2.gff

Will create file1.no-missing-scores.gff and file2.no-missing-scores.gff

=cut


