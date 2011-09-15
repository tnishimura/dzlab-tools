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
use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "output|o=s" => \my $output,
    "help|h" => \my $help,
);
pod2usage(-verbose => 99) if ($help);  

if ($output){
    open my $fh, '>', $output;
    select $fh;
}

my $parser = GFF::Parser->new(file => \*ARGV);
while (defined(my $gff = $parser->next())){
    my $s = $gff->score();
    my $start = $gff->start();
    my $end = $gff->end();
    if (defined $s && defined $start && defined $end){
        $gff->score($s / ($end - $start + 1));
    }
    say $gff->to_string;
}

=head1 NAME

gff_average.pl - replace scores with score/(end-start+1) to get per-base scores

=head1 SYNOPSIS

Usage examples:

 gff_average.pl -o output.gff input.gff
 gff_average.pl input.gff > output.gff
 some_other_script.pl | gff_average.pl > output.gff

=head1 OPTIONS

=over

=item --help | -h

=item  -o <file> | --output <file>

=back

=cut
