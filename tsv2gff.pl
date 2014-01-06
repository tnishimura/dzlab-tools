#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use IO::File;
use IO::Handle;
use Pod::Usage;
use Getopt::Long;
use Text::CSV_XS;

my $result = GetOptions (
    "sequence-column|c=i" => \(my $sequence_column),
    "start-column|s=i"    => \(my $start_column),
    "end-column|e=i"      => \(my $end_column),
    "ignore-column|g=i"   => \(my $ignore_column),
    "width|w=i"           => \(my $width),
    "sep-char=s"          => \(my $sep_char = "\t"),
    "input|i=s"           => \(my $input_file),
    "zero-base|z"         => \(my $zero_base=1),
    "no-skip|k"           => \(my $no_skip),
);

pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || ! defined $sequence_column || ! defined $start_column || ! $input_file);  

if (! ($width xor $end_column)){
    die "you must specify either --width/-w OR --end-column/-e, not both.";
}

if ($zero_base){
    $sequence_column--;
    $start_column--;
    $end_column-- if defined $end_column;
}

open my $fh, '<:crlf', $input_file;
my $csv = Text::CSV_XS->new({ eol => "\n", sep_char => $sep_char, });
my %handles;

#######################################################################
# header

my $header = $csv->getline($fh);

for my $i (0 .. $#{$header}) {
    if (is_value_column($i)){
        my $colname = $header->[$i];
        if ($colname !~ /[[:ascii:]]+/){
            die "$colname bad";
        }
        my $fn = make_file_name($input_file, $colname);
        if (exists $handles{$i}){
            die "multiple columns named $colname exists?";
        }

        $handles{$i} = IO::File->new($fn, 'w');
    }
}

#######################################################################
# body

while (my $row = $csv->getline ($fh)) {
    my $seq = $row->[$sequence_column];
    my $start = $row->[$start_column];
    my $end = defined $width ? $start + $width : $row->[$end_column];

    for my $i (0 .. $#{$header}) {
        if (is_value_column($i)){
            my $score = $row->[$i];
            if (defined $score and $score ne '' and $score ne '.'){
                $handles{$i}->print(
                    join("\t", $seq, '.', '.', $start, $end, $row->[$i], '.', '.', '.',) . "\n"
                );
            }
            elsif ($no_skip){
                $handles{$i}->print(
                    join("\t", $seq, '.', '.', $start, $end, '.', '.', '.', '.',) . "\n"
                );
            }
        }
    }
    
}

#######################################################################
# finalize

$csv->eof or $csv->error_diag ();
close $fh;
close $_ for values %handles;

#######################################################################

sub make_file_name{
    my ($input_file, $header_name) = @_;
    $header_name =~ s/[^\w]//g;
    if ($input_file =~ s/\.(\w+)$//){
        return "$input_file.$header_name.gff";
    }
    else{
        return "$input_file.$header_name.gff";
    }
}

sub is_value_column{
    my $col = shift;
    return (
        $col != $sequence_column
            and 
        $col != $start_column 
            and 
        (! defined $end_column or $col != $end_column) 
            and 
        (! defined $ignore_column or $col != $ignore_column)
    );
}

=head1 NAME

tsv2gff.pl - split a mashup file/delimiter-separated file (CSV, TSV) into gff files.

=head1 SYNOPSIS

Split a mashup file (col #1 is seqid, #2 is start position, each line is 50bp) into gff files:

 tsv2gff.pl -i input.mashup -c 1 -s 2 -w 50 -i input.txt

Split a tab-separated file where col #1 is seqid, #2 is start position, #3 is end:

 tsv2gff.pl -i input.mashup -c 1 -s 2 -e 3 -i input.txt

=head1 OPTIONS

=over

=item -i filename | --input filename

Input file, required.

=item -c # | --sequence-column #

Column number of the sequence identifier (chr1, chr2, etc). Required.

=item -s # | --start-column #

Column number of start position column. Required.

=item -e # | --end-column #

Column number of end position column.  You can specify this OR --end-column. 

=item -w # | --width #

Width of each row.  You can specify this OR --end-column. 

=item -g # | --ignore-column #

Column number of column to ignore.

=item --sep-char ,

Separation character, defaults to TAB.

=item -k | --no-skip

Do not skip values which are empty.

=item -z | --zero-base

Columns are base 0, instead of base 1 (the default).

=back

=cut

