#!/usr/bin/env perl
package main;
use strict;
use warnings FATAL => "all";
use 5.010_000;
use FindBin;
use lib "$FindBin::Bin/lib";
use autodie;
use Data::Dumper;
use File::Temp qw/tempdir tempfile/;
use File::Basename;
use Pod::Usage;
use Log::Log4perl qw/:easy/;
use Launch qw/cast/;
use List::Util qw/max/;
use Pod::Usage;
use Getopt::Long;
use Cwd qw/getcwd/;
use File::Basename qw/basename dirname/;
use File::Path qw/make_path remove_tree/;
use File::Spec::Functions qw/rel2abs canonpath catdir catfile updir/;
use File::Copy;


END {close STDOUT}
$| = 1;



#######################################################################
# utilities

sub cleanse{
    my $name = shift;
    $name =~ s/^\s+//;
    $name =~ s/\s+$//;
    return $name;
}

# parse_key_specs("1,2,3,5n")
#  => [[col#, 't' | 'n'], ...]
sub parse_key_specs{
    my $colspec = shift;

    my @split = map { cleanse $_ } split /,/, $colspec;

    for (@split){ 
        if (! /^\d+n?$/){
            die "'$_' is not a valid key specification"
        }
    }

    return [map {
        my $col = $_;
        $col =~ s/n$// && [$col, 'n'] || [$col, 't'];
    } @split];
}

sub parse_val{
    my $colspec = shift;

    my @split = map { cleanse $_ } split /,/, $colspec;

    for (@split){ 
        if (! /^\d+$/){
            die "$_ is not a valid value specification"
        }
    }

    return @split;
}

sub extract_columns{
    my ($parts, $columns) = @_;
    return map {
        $parts->[$_ - 1];
    } @$columns;
}

sub sort_file_in_place{
    my ($file, $colspec) = @_;
    my (undef, $tempfile) = tempfile("$file.XXXXX");
    cast("sort", "-S", "5%", "-f", 
        (
            map {
                my $colnum = $_->[0];
                my $coltype  = $_->[1];
                my $spec = "-k$colnum,$colnum";
                if ($coltype eq 'n') {
                    $spec .= "n";
                }
                elsif ($coltype ne 't'){
                    die "";
                }
                $spec
            } @$colspec
        ),
        $file,
        '-o',
        $tempfile);

    rename $tempfile, $file;
}

# parse_nicknames ["file1.txt", "file2.txt", ...], "nick1,nick2,..."
# => [["file1.txt" => 'nick1'], ["file2.txt" => 'nick2', ...]], 
sub parse_nicknames{
    my ($filenames, $nicknames) = @_;

    my @nicks = defined $nicknames ? (map { cleanse $_ } split ',', $nicknames) : ();

    if (@nicks > 0){
        if (@nicks != @$filenames){
            die "If you're going to specific nicks, number of nicks must match number of input files";
        }
        else{
            return [map { [$filenames->[$_], $nicks[$_]] } (0 .. $#nicks)];
        }
    }
    else{
        return [map { [$_, basename($_)] } @$filenames];
    }
}

sub make_comparator{
    my $keycolspec = shift;
    my $maxcol = max map { $_->[0] } @$keycolspec;

    return sub{
        my ($left_parts, $right_parts) = @_;
        my $left_size  = @$left_parts;
        my $right_size = @$right_parts;

        die "uneven parts size" if $left_size != $right_size;

        die "not enough columns" if ($maxcol > $left_size || $maxcol > $right_size);

        #say Dumper "make_comparator", $left_parts, $right_parts;

        for my $pair (@$keycolspec) {
            my ($colnum, $type) = @$pair;
            my $left = $left_parts->[$colnum - 1];
            my $right = $right_parts->[$colnum - 1];
            my $cmp;
            if ($type eq 't'){
                $cmp = $left cmp $right;
            }
            elsif ($type eq 'n'){
                $cmp = $left <=> $right;
            }

            if ($cmp != 0){
                return $cmp;
            }
        }
        return 0;

    }
}

{
    my $sep = "\t";
    sub make_split_iterator{
        my $filename = shift;
        my $numcol;
        my $numread = 0;
        my @rewind;
        open my $fh, '<', $filename;
        return sub{
            my $cmd = shift;
            if ($cmd){
                if ($cmd eq 'eof'){
                    return eof $fh;
                }
                elsif ($cmd eq 'numread'){
                    return $numread;
                }
                elsif ($cmd eq 'rewind'){
                    push @rewind, @_;
                }
            }
            elsif (! eof $fh){
                if (@rewind){
                    return pop @rewind;
                }
                else{
                    my $line = readline($fh);
                    $line =~ tr/\n\r//d;
                    my @split = split /$sep/, $line, -1;
                    $numcol //= scalar @split;
                    if ($numcol != scalar @split){
                        die "uneven column counts in $filename";
                    }
                    $numread++;
                    return \@split;
                }
            }
            else{
                return;
            }
        };
    }
}

sub ndots{
    my $n = shift;
    return join "\t", ('.') x $n;
}

sub run{
    my $result = GetOptions (
        "output|o=s"    => \my $opt_output,
        "nicknames|n=s" => \my $opt_nicknames,
        "values|v=s"   => \my $opt_columns,
        "keys|k=s"      => \my $opt_keys,
    );
    if (!$result || ! $opt_output || ! $opt_columns || ! $opt_keys){
        say "usage: generic_mashup.pl ...";
        exit 1;
    }

    my @opt_input = @ARGV;

    my $key_specs       = parse_key_specs($opt_keys);
    my @key_columns     = map { $_->[0] } @$key_specs;
    my @val_columns     = parse_val($opt_columns); 
    my $comparator      = make_comparator($key_specs);
    my $files_and_nicks = parse_nicknames(\@opt_input, $opt_nicknames);

    say Dumper \@opt_input, $files_and_nicks;

    for my $pair (@$files_and_nicks) {
        my ($input_file, $nick) = @$pair;

        sort_file_in_place($input_file, $key_specs);

        # output to tempfile
        my (undef, $tempout_file) = tempfile("$opt_output.XXXXX");

        open my $tempout, '>', $tempout_file;
        my $input_parser = make_split_iterator($input_file);

        if (-f $opt_output){
            my $mashup_parser = make_split_iterator($opt_output);
            my $num_mashup_vals;
            while (defined(my $mashup_line = $mashup_parser->())){
                my @mashup_keys = extract_columns($mashup_line, \@key_columns);
                my @mashup_vals = extract_columns($mashup_line, \@val_columns);
                $num_mashup_vals //= @mashup_vals;

                if ($mashup_parser->('numread') == 1){
                    say $tempout join "\t", @mashup_keys, @mashup_vals, map { "${nick}_$_" } @val_columns;
                }
                else{
                    if (! $input_parser->('eof')){
                        INPUT:
                        while (defined(my $input_line = $input_parser->())){
                            my @input_keys = extract_columns($input_line, \@key_columns);
                            my @input_vals = extract_columns($input_line, \@val_columns);
                            

                            my $cmp = $comparator->($input_line, $mashup_line);

                            # input behind mashup
                            if ($cmp == -1){
                                say $tempout join "\t", @input_keys, ndots($num_mashup_vals), @input_vals;
                            }
                            # input matches mashup
                            elsif ($cmp == 0){
                                say $tempout join "\t", @mashup_keys, @mashup_vals, @input_vals;
                                last INPUT; # and get next mashup line
                            }
                            # input ahead of mashup
                            elsif ($cmp == 1){
                                say $tempout join "\t", @mashup_keys, @mashup_vals, ndots(scalar(@input_vals));
                                $input_parser->('rewind', $input_line);
                                last INPUT; # and get next mashup line
                            }
                        }
                    }
                    # excess lines in mashup file
                    else{
                        say $tempout join "\t", @mashup_keys, @mashup_vals, ndots(scalar(@val_columns));
                    }
                }
            }

            # excess lines in input file
            while (defined(my $input_line = $input_parser->())){
                my @input_keys = extract_columns($input_line, \@key_columns);
                my @input_vals = extract_columns($input_line, \@val_columns);

                say $tempout join "\t", @input_keys, ndots($num_mashup_vals), @input_vals;
            }
        }

        # no existing output file
        else{
            say $tempout join "\t", (map { "key_$_" } @key_columns), (map { "${nick}_$_" } @val_columns);
            while (defined(my $parts = $input_parser->())){
                say $tempout join "\t", extract_columns($parts, \@key_columns), extract_columns($parts, \@val_columns);
            }
        }
        close $tempout;
        rename $tempout_file, $opt_output;
    }
}

run() unless caller();

# for my $pair (@files) {
#     my ($input_file, $nick) = @$pair;
# 
#     # input
#     my $parser = GFF::Parser->new(file => $input_file, normalize => 0);
# 
#     # existinng output. append.
#     if (-f $opt_output){
#         open my $output_read, '<', $opt_output;
#         my $numcells;
#         while (defined(my $line = <$output_read>)){
#             chomp $line;
#             # -1 to split b/c when stata mangles output file, dots are turned to blank...
#             my ($seq, $coord, @cells) = split /\t/, $line, -1; 
#             $numcells //= @cells;
# 
#             # header
#             if ($seq =~ /^Sequence/i){
#                 say $tempout join "\t", $seq, $coord, @cells, map { $nick . "_" . $_ } @columns;
#             }
#             else {
#                 if (! $parser->done()){
#                     PARSER:
#                     while (defined (my $gff = $parser->next())){
#                         # gff behind output. new line with gff.
#                         if (gff_lessthan($gff->sequence, $gff->start, $seq, $coord)){
#                             say $tempout join "\t", $gff->sequence, $gff->start, (map { '.' } @cells), get_columns($gff, @columns);
#                         }
#                         # gff past output. new line with output. and move to next output
#                         elsif (gff_greaterthan($gff->sequence, $gff->start, $seq, $coord)){
#                             say $tempout join "\t", $seq, $coord, @cells, map { '.' } @columns;
#                             $parser->rewind($gff);
#                             last PARSER;
#                         }
#                         # coincides
#                         elsif (gff_equal($gff->sequence, $gff->start, $seq, $coord)){
#                             say $tempout join "\t", $seq, $coord, @cells, get_columns($gff, @columns);
#                             last PARSER;
#                         }
#                     }
#                 }
#                 else {
#                     # print excess outputs on new line
#                     say $tempout join "\t", $seq, $coord, @cells, map { '.' } @columns;
#                 }
#             }
#         }
#         # get any remaining that are past the end of output:
#         while (defined (my $gff = $parser->next())){
#             say $tempout join "\t", $gff->sequence, $gff->start, (map { '.' } (1 .. $numcells)), get_columns($gff, @columns);
#         }
# 
#         close $output_read;
#     }
# 
#     # new output file
#     else {
#         say $tempout join "\t", "Sequence", "Coord", map { $nick . "_" . $_ } @columns;
#         while (defined(my $gff = $parser->next())){
#             say $tempout join "\t", $gff->sequence, $gff->start, get_columns($gff, @columns);
#         }
#     }
#     close $tempout;
# 
#     rename $tempout_file, $opt_output;
# }
