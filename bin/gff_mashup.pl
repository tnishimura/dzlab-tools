#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use FindBin;
use lib "$FindBin::Bin/../lib";
use autodie;
use Data::Dumper;
use File::Temp qw/tempdir tempfile/;
use File::Basename;
use Getopt::Euclid qw( :vars<opt_> );
use GFF::Parser;
use Pod::Usage;
use Log::Log4perl qw/:easy/;
use Launch;

END {close STDOUT}
$| = 1;

Log::Log4perl->easy_init({ level => $DEBUG, layout => '%d{HH:mm:ss} %.1p > %m%n' });
my $logger = get_logger();

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || ! @opt_input || ! $opt_output;

my @columns = split ',', $opt_columns;
my @nicks = defined $opt_nicknames ? split ',', $opt_nicknames : ();

if (@nicks > 0 && @nicks != @opt_input){
    die "If you're going to specific nicks, number of nicks must match number of input files";
}

my @files = scalar @nicks ? 
(map { [$opt_input[$_], $nicks[$_]] }   (0 .. $#opt_input)) :
(map { [$opt_input[$_], basename($opt_input[$_])] } (0 .. $#opt_input));

for my $pair (@files) {
    my ($input_file, $nick) = @$pair;

    # sort unless --no-sort
    unless ($opt_no_sort){
        my (undef, $tempfile) = tempfile("$input_file.XXXXX");

        $logger->info("Sorting $input_file to $tempfile");
        launch(qq/sort -f -k1,1 -k4,4n "$input_file" -o "$tempfile"/);

        $logger->info("Overwriting $input_file with sorted version. (don't quit script right now)");
        rename($tempfile, $input_file);

        $logger->info("Done");
    }
    
    # output to tempfile
    my (undef, $tempout_file) = tempfile("$opt_output.XXXXX");
    open my $tempout, '>', $tempout_file;

    # input
    my $parser = GFF::Parser->new(file => $input_file, normalize => 0);

    # existinng output. append.
    if (-f $opt_output){
        open my $output_read, '<', $opt_output;
        my $numcells;
        while (defined(my $line = <$output_read>)){
            chomp $line;
            # warn "output line: $line";
            # -1 to split b/c when stata mangles output file, dots are turned to blank...
            my ($seq, $coord, @cells) = split /\t/, $line, -1; 
            $numcells //= @cells;
            my $output_line_cells_printed = 0;

            # header
            if ($seq =~ /^Sequence/i){
                say $tempout join "\t", $seq, $coord, @cells, map { $nick . "_" . $_ } @columns;
            }
            else {
                if (! $parser->done()){
                    PARSER:
                    while (defined (my $gff = $parser->next())){
                        # gff behind output. new line with gff.
                        if (gff_lessthan($gff->sequence, $gff->start, $seq, $coord)){
                            # warn "<";
                            say $tempout join "\t", $gff->sequence, $gff->start, (map { '.' } @cells), get_columns($gff, @columns);
                        }
                        # gff past output. new line with output. and move to next output
                        elsif (gff_greaterthan($gff->sequence, $gff->start, $seq, $coord)){
                            # warn ">";
                            say $tempout join "\t", $seq, $coord, @cells, map { '.' } @columns;
                            $output_line_cells_printed = 1;
                            $parser->rewind($gff);
                            last PARSER;
                        }
                        # coincides
                        elsif (gff_equal($gff->sequence, $gff->start, $seq, $coord)){
                            # warn "=";
                            say $tempout join "\t", $seq, $coord, @cells, get_columns($gff, @columns);
                            $output_line_cells_printed = 1;
                            last PARSER;
                        }
                    }
                }
                if (! $output_line_cells_printed){
                    say $tempout join "\t", $seq, $coord, @cells, map { '.' } @columns;
                }
            }
        }
        # get any remaining that are past the end of output:
        while (defined (my $gff = $parser->next())){
            say $tempout join "\t", $gff->sequence, $gff->start, (map { '.' } (1 .. $numcells)), get_columns($gff, @columns);
        }

        close $output_read;
    }

    # new output file
    else {
        say $tempout join "\t", "Sequence", "Coord", map { $nick . "_" . $_ } @columns;
        while (defined(my $gff = $parser->next())){
            say $tempout join "\t", $gff->sequence, $gff->start, get_columns($gff, @columns);
        }
    }
    close $tempout;

    rename $tempout_file, $opt_output;
}

sub get_columns{
    my ($gff, @columns) = @_;
    my @accum;
    my ($c, $t);

    for my $col (@columns) {
        given ($col){
            when ('c+t'){
                $c = $gff->get_column('c');
                $t = $gff->get_column('t');
                if (! defined $c && ! defined $t){
                    push @accum, '.';
                }
                else{
                    push @accum, ($c//0) + ($t//0);
                }
            }
            when ('methyl'){
                $c = $gff->get_column('c');
                $t = $gff->get_column('t');
                if (! defined $c && ! defined $t){
                    push @accum, '.';
                }
                else{
                    my $total = ($c//0) + ($t//0);
                    push @accum, $total > 0 ? sprintf("%.4f", ($c//0) / $total) : 0;
                }
            }
            default{
                push @accum, $gff->get_column($_) // '.';
            }
        }
    }
    return @accum;
}

sub gff_lessthan{
    my ($seq1, $coord1, $seq2, $coord2) = @_;
    my $cmp = lc $seq1 cmp lc $seq2  || $coord1 <=> $coord2;
    return $cmp == -1;
}
sub gff_greaterthan{
    my ($seq1, $coord1, $seq2, $coord2) = @_;
    my $cmp = lc $seq1 cmp lc $seq2  || $coord1 <=> $coord2;
    return $cmp == 1;
}
sub gff_equal{
    my ($seq1, $coord1, $seq2, $coord2) = @_;
    return lc $seq1 eq lc $seq2  && $coord1 == $coord2;
}

=head1 NAME

gff_mashup.pl - Extract columns from GFF files and create a tab delimited file with columns from each file.

=head1 SYNOPSIS

Usage examples:

Extract the single-c scores from input gff files (not that 'c+t,methyl' is
default for -c):

 gff_mashup.pl -o output.txt input1.gff input2.gff input3.gff 

Extract C's, T's, and strand from inputs by specifying a custom column
specification:

 gff_mashup.pl -c 'c,t,strand' -o output.txt input1.gff input2.gff input3.gff 

Extract C's, T's, and strand from inputs, nickname the columns WT, Mut1, Mut2
by specifying custom nicknames:

 gff_mashup.pl -n 'WT,Mut1,Mut2' -c 'c,t,strand' -o output.txt input1.gff input2.gff input3.gff 

NOTE: You can specify an existing output file and the script will APPEND rather than OVERWRITE
the scores appropriately.  In other words, the commands:

 gff_mashup.pl -o output.txt input1.gff input2.gff 
 gff_mashup.pl -o output.txt input3.gff 

and 

 gff_mashup.pl -o output.txt input1.gff input2.gff input3.gff 

are equivalent.

=head1 OPTIONS

=over

=item  -o <out> | --output <out>

=item  -c <columns> | --columns <columns>

Comma separate list of columns you want to extract from input gff files.  Can
be any of the standard nine columns: sequence, source, feature, start, end,
score, strand, frame, or attributes.  OR It can be any key in the attributes
field (so if you wanted to extract "AT12345" from "ID=AT12343", you can specify
ID as a column.  OR you can specify 'c+t' or 'methyl', which calculates the c+t
or methyl scores.

=for Euclid
    columns.default:     'c+t,methyl'

=item  <input>...

Input files.

=for Euclid
    input.type:        readable

=item -n <nicknames> | --nicknames <nicknames>

Comma separated list of nicknames for files.  Used as column headers. Number of
nicknames must match number of input files.

=item  -k  | --no-sort

Assuming input gff files are sorted.

=item --help | -h

=back

=cut
