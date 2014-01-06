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
use Text::CSV_XS;
use Scalar::Util qw/looks_like_number/;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || !$opt_dictionary || (! $opt_input && -t STDIN);


##

my @value_columns = $opt_columns ? map { 
                                       die "$opt_columns not csv of numbers?" if ! looks_like_number($_); 
                                       $_ - 1 
                                   } split /,/, $opt_columns
                                 : (4,5);
my $sequence_column = $opt_seq_column ? $opt_seq_column - 1 : 0;
my $is_gff = ! defined $opt_columns && ! defined $opt_seq_column;


my $input_fh;
if (ref $opt_input eq 'GLOB'){
    $input_fh = $opt_input;
}
else{
    open $input_fh, '<:crlf', $opt_input ;
}

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}


my $config = LoadFile("$FindBin::Bin/conf/$opt_dictionary.alignment");
my $dictionary = $config->{alignment};
#die Dumper $dictionary;

sub convert{
    my ($dictionary, $seq, $coord) = @_;
    $seq = uc $seq;
    my $alignment = $dictionary->{$seq};
    #say Dumper $alignment;
    if (defined $alignment){
        # linear search... ok for now.
        for my $island (@$alignment) {
            if ($island->[0] <= $coord && $coord <= $island->[1]){
                return $island->[2] + ($coord - $island->[0]);
            }
        }
        warn("no appropriate island found for $seq, $coord, line $.") if $opt_debug;
        return;
    }
    else {
        warn("$seq not found in dictionary") if $opt_debug;
        return $coord;
    }
}


if ($is_gff){
    my $parser = GFF::Parser->new(file => $input_fh, normalize => 0);

    while (defined(my $gff = $parser->next())){
        # warn($gff->to_string) if $opt_debug;
        my $start = convert $dictionary, uc($gff->sequence()), $gff->start();
        my $end   = convert $dictionary, uc($gff->sequence()), $gff->end();
        if ($start && $end){
            $gff->start($start);
            $gff->end($end);
            say $gff->to_string;
        } 
        # else {
        #     warn("could not convert: " . $gff->to_string) if $opt_debug;
        #     # say $gff->to_string;
        # }
    }
}
else{
    my $csv = Text::CSV_XS->new({
        eol      => "\n",
        sep_char => "\t",
    });
    
    ROWLOOP:
    while (my $row = $csv->getline($input_fh)) {
        my $seqid = $row->[$sequence_column];
        for my $vcol (@value_columns) {
            if (! looks_like_number($row->[$vcol])){
                increment_error("column %d ('%s') of line '%s' is not numeric? skipping...", $vcol + 1, $row->[$vcol], join "\t", @$row);
                next ROWLOOP;
            }
            $row->[$vcol] = convert $dictionary, uc($seqid), $row->[$vcol];
            if (! defined $row->[$vcol]){
                next ROWLOOP;
            }
        }
        say join "\t", @$row;
    }
    $csv->eof or $csv->error_diag ();
}

if (ref $opt_input ne 'GLOB'){
    close $input_fh;
}

{
    my $error = 0;
    sub increment_error{
        my ($format, @args) = @_;
        if (++$error > $opt_max_error){
            die(sprintf("$format\ndying b/c too many\n",@args));
        }
        else {
            warn(sprintf("$format\n",@args));
        }
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

=for Euclid
    file.default:     '-'

=item  <input> 

Input file.

=for Euclid
    input.type:        readable

=item  -d <dict> | --dictionary <dict>

Conversion dictionary.  These are files found in the coord/ directory in the
install path (c:\dzlab-tools\coord on windows).  The '.alignment' prefix is
optional. Currently available are:

 rice-5.0-to-6.1
 rice-6.1-to-5.0
 at-tair8-to-tair10
 at-tair10-to-tair8

=item  -s <col> | --seq-column <col>

Column position of sequence name.  Default to 1 for the seqid column in a GFF.  

=for Euclid
    col.type:        readable
    col.type:        int, col >= 1 
    col.type.error:  <col> must be greater than 1. 

=item  -c <col> | --columns <col>

Comma separated list of columns of convert.  Default is '4,5' for the start and
end columns in a GFF.  Note that there cannot be a space in the list, so '4, 5'
or '4 ,5' are invalid and will error. 

=item  -e <error_count> | --max-error <error_count>

=for Euclid
    error_count.default:     10
    error_count.type:        int, error_count >= 1 
    error_count.type.error:  <error_count> must be greater than 1 

=item  --debug

=item -h | --help

=back

=cut

