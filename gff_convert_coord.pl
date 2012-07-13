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
use GFF::Parser;
use Log::Log4perl qw/:easy/;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || ! $opt_output || !$opt_input || !$opt_dictionary;

Log::Log4perl->easy_init( { 
    level    => $TRACE,
    layout   => '%p> (%L) %M - %m%n',
    (defined $opt_logfile ? (file => ">$opt_logfile") : ()),
} );

my $logger = get_logger();

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}

my $config = LoadFile("$FindBin::Bin/conf/$opt_dictionary.alignment");
#die "$FindBin::Bin/conf/$opt_dictionary.alignment";
my $dictionary = $config->{alignment};
#die Dumper $dictionary;

# $dictionary == {seq => [[from_start,from_end, to_start, to_end] ... ]}

sub convert{
    my ($dictionary, $seq, $coord) = @_;
    $seq = uc $seq;
    my $alignment = $dictionary->{$seq};
    #say Dumper $alignment;
    if (defined $alignment){
        #$logger->trace("found alignment");
        # linear search... ok for now.
        for my $island (@$alignment) {
            if ($island->[0] <= $coord && $coord <= $island->[1]){
                return $island->[2] + ($coord - $island->[0]);
            }
        }
        $logger->debug("no appropriate island found for $seq, $coord");
    }
    else {
        # $logger->debug("$seq not found in dictionary");
        return $coord;
    }
    return;
}

#say Dumper $dictionary;
#say convert $dictionary, "CHR03", 1000;

my $parser = GFF::Parser->new(file => $opt_input, normalize => 0);

while (defined(my $gff = $parser->next())){
    #$logger->trace($gff->to_string);
    my $start = convert $dictionary, uc($gff->sequence()), $gff->start();
    my $end   = convert $dictionary, uc($gff->sequence()), $gff->end();
    if ($start && $end){
        $gff->start($start);
        $gff->end($end);
        say $gff->to_string;
    } else {
        $logger->info("could not convert: " . $gff->to_string);
        # say $gff->to_string;
    }
}

=head1 NAME
 
gff_convert_coord.pl - convert coordinates according to a dictionary.
 
=head1 SYNOPSIS

Convert column 4 and 5 in input.gff to output.gff according to rice-5.0-to-6.1:

 gff_convert_coord.pl -d rice-5.0-to-6.1 -c 4,5 -o output.gff input.gff

=head1 OPTIONS

=over

=item -o <file> | --output <file>

Output file.

=item  <input> 

Input file.

=for Euclid
    input.type:        readable

=item  -l <log> | --logfile <log>

Put error warnings in this file.  Default to the screen.

=item  -d <dict> | --dictionary <dict>

Conversion dictionary.  These are files found in the coord/ directory in the
install path (c:\dzlab-tools\coord on windows).  The '.alignment' prefix is
optional. Currently available are:

 rice-5.0-to-6.1
 rice-6.1-to-5.0
 at-tair8-to-tair10
 at-tair10-to-tair8

=item  -c <col> | --columns <col>

Comma separated list of columns of convert.  Default is '4,5' for the start and
end columns in a GFF.  Note that there cannot be a space in the list, so '4, 5'
or '4 ,5' are invalid and will error. 

=for Euclid
    col.default:     '4,5'

=item -h | --help

=back

=cut

