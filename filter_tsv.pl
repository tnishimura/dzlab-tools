#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use lib "$FindBin::Bin/Tree-Range/lib";
use GFF::Tree;
use GFF::Parser;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use File::Temp qw/tempfile/;
use Scalar::Util qw/looks_like_number/;

END {close STDOUT}
$| = 1;


pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if ($opt_in_place && $opt_output) 
|| ! $opt_start_column
|| ($opt_end_column && $opt_window_size)
|| ! $opt_annotation 
|| ! $opt_input 
|| $opt_help;


#######################################################################
# Setup input and output. maybe this can be modularized?

my $input_fh;
my $output_fh;
if ($opt_in_place){
    my (undef, $tempfile) = tempfile("$opt_input.XXXXX", UNLINK => 0);
    rename $opt_input, $tempfile;
    open $input_fh, '<', $tempfile;
    open $output_fh, '>', $opt_input;
}
elsif ($opt_output){
    open $input_fh, '<', $opt_input;
    open $output_fh, '>', $opt_output;
}
else{
    open $input_fh, '<', $opt_input;
    $output_fh = \*STDOUT;
}

#######################################################################
# Header?

my $nick = $opt_nickname // basename($opt_input);
if ($opt_skip_header){
    my $header = <$input_fh>;
    chomp $header;
    say $output_fh "$header\t$nick";
}

my $indicator = $opt_nick_in_cell ? $nick : 1;

#######################################################################
# Loop

my $tree = GFF::Tree->new(file => $opt_annotation, normalize => 1);

while (defined(my $line = <$input_fh>)){
    chomp $line;
    my @F = split /\t/, $line;
    my $seq = uc $F[$opt_seq - 1];
    my $start = $F[$opt_start_column - 1];
    my $end = $opt_end_column ? $F[$opt_end_column - 1] : 
              $opt_window_size ? ($start + $opt_window_size) : 
              $start;
    if (defined $seq && looks_like_number($start) && looks_like_number($end)){
        my $overlap = $tree->search_overlap($seq, $start, $end);
        say $output_fh join "\t", @F, ($overlap ? $indicator : 0);
    }
    else {
        say $output_fh join "\t", @F, '.';
    }
}


=head1 NAME

filter_tsv.pl - compare start/end coordinates of a tab separated file (tsv,
including GFF's) and add a column indicating whether entry is in annotation
file.

=head1 SYNOPSIS

Filtering a GFF file, in place:

 filter_tsv.pl -s 4 -e 5 -i -a genes.gff input.gff

Filtering a gff_mashup.pl output file, in place:

 filter_tsv.pl -s 2 -i -h -a genes.gff -n genes input.gff

Filtering a gff_mashup.pl output file, outputing to output.txt:

 filter_tsv.pl -s 2 -i -h -a genes.gff -n genes -o output.txt input.gff


=head1 OPTIONS

=over

=item  -seq <seq_column> 

Sequence name column.  Should be 1 (the default) for both GFF and GFF mashup
files.

=for Euclid
    seq_column.default:     1

=item  -s <coord> | --start-column <coord>

Column position of the start column (with the leftmost column being column 1,
not 0).  For GFF file, use 4. For GFF mashup output, use 2.

=item  -e <coord> | --end-column <coord>

Column position of the end column.  For GFF file, use 5. For GFF mashup, don't use this option, since 
it won't have an end column.  Use -w instead.  -e and -w are mutually exclusive.

=item  -w <size> | --window-size <size>

If specified instead of -e, assume that window ends are of "start_coord +
window_size - 1". This is useful for GFF mashup output which only have a start
column, or when you want to force a single-c GFF file to be treated as a window
X file.  -e and -w are mutually exclusive.

=item  -a <gff_annotation_file> | --annotation <gff_annotation_file>

GFF annotation file.

=for Euclid
    gff_annotation_file.type:        readable

=item  -n <nickname> | --nickname <nickname>

Optional name for new column's header.  If omitted, the file name is used.

=item  -i | --in-place

=item  -o <file> | --output <file>

=item  -h | --skip-header 

Skip a single line at the beginning of the input file.  For files with headers,
like GFF mashup output.

=item -c | --nick-in-cell

If this option is used, the cells in the new column will be the nickname, not '1'.

=item  <input> 

=for Euclid
    input.type:        readable

=item --help | -h

=back

=cut
