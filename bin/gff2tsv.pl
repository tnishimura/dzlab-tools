#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/all/;
use List::Util qw/max min/;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;
use Pod::Usage;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../lib";
#use GFF::Parser::Attributes;
use GFF::Parser;
use Log::Log4perl qw/:easy/;
Log::Log4perl->easy_init( { 
    level    => $DEBUG,
    #file     => ">run.log",
    layout   => '%d{HH:mm:ss} %p> (%L) %M - %m%n',
} );
my $logger = get_logger();

my $output = q{-};
my $help;
my $result = GetOptions (
    "output|o=s"    => \$output,
    "help"          => \$help,
);
pod2usage(-verbose => 1) if (!$result || $help || ! @ARGV);

unless ($output eq '-'){
    close STDOUT;
    open STDOUT, '>', $output or die "can't open $output for writing";
}

my @mains = qw/sequence source feature start end score   strand frame   attribute_string/;
my %seen;

$logger->info("Finding out columns...");

# read through all files once, recording all attribute names
foreach my $file (@ARGV) {
    my $parser = GFF::Parser->new(file => $file);
    while (my $gff = $parser->next()){
        my @attributes = $gff->list_attribute();
        @seen{@attributes} = map { 1 } @attributes;
    }
}

my @columns = (qw/sequence source feature start end score strand frame/, sort keys %seen);

# read through them again, this time printing.

$logger->info("Converting...");

say join "\t",@columns;
foreach my $file (@ARGV) {
    my $parser = GFF::Parser->new(file => $file);
    while (my $gff = $parser->next()){
        say join "\t", map { $gff->get_column($_) // '.' } @columns;
    }
}

=head1 NAME

gff2tsv.pl - convert a GFF file to TSV, tab seperated file

=head1 SYNOPSIS

 gff2tsv.pl -o output-tsv.txt input1.gff input2.gff ...

=head1 DESCRIPTION

The first 8 columns are always seqname, source, feature, start, end, score, strand, frame.  The
ninth and further columns are the attributes split up into individual column.  First line of output
is the column names.

=cut

