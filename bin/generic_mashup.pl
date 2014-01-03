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
use Pod::Usage;
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

# extract_columns([qw/col1 col2 col3 .../], [3,2]) => qw/col3 col2/
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
    my $numkeys_expected = @$keycolspec;
    my $maxcol = max map { $_->[0] } @$keycolspec;

    return sub{
        my ($input_keys, $mashup_keys) = @_;

        die "bug: unexpected number of mashup_keys. report to programmer" if @$mashup_keys != $numkeys_expected;
        die "bug: unexpected number of input_keys. report to programmer" if @$input_keys != $numkeys_expected;

        #say Dumper "make_comparator", $input_line, $mashup_keys;

        for my $i (0 .. $#$keycolspec) {
            #for my $pair (@$keycolspec) {
            my $colnum = $keycolspec->[$i][0];
            my $type = $keycolspec->[$i][1];

            my $left = $input_keys->[$i];
            my $right = $mashup_keys->[$i];

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

# line-splitting file iterator.
{
    my $sep = "\t";
    sub make_split_iterator{
        my $filename = shift;
        my $numcol;
        my $numread = 0;
        my @rewind; # fifo of put-backs
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

# create list of n dots
sub ndots{
    my $n = shift;
    return join "\t", ('.') x $n;
}

sub run{
    $0 = basename($0);
    my $result = GetOptions (
        "output|o=s"    => \my $opt_output,
        "nicknames|n=s" => \my $opt_nicknames,
        "values|v=s"    => \my $opt_columns,
        "keys|k=s"      => \my $opt_keys,
    );
    if (!$result || ! $opt_output || ! $opt_columns || ! $opt_keys){
        say <<"END";
$0 - mashup tab-separated files using certain columns as keys and other 
     columns as values

     --output    -o  Output file.
     --keys      -k  Comma separated list of columns to use as key.  
                     Append numerical columns with a 'n'.
     --values    -v  Comma separated list of columns to use as values.  
                     Append numerical columns with a 'n'.
     --nicknames -n  Comma separated list of nicknames of files. Optional.

usage:

use column 1 (as text) and 2 (as numerical) as keys to mashup column 3:

  $0 -o output.txt --keys '1,2n' --values '3' input1.txt input2.txt input3.txt

use column 1 (as text) as key to mashup column 3, 4:

  $0 -o output.txt --keys '1' --values '3,4' input1.txt input2.txt input3.txt

use column 3, then 2 keys to mashup column 3, with file nicknames:

  $0 -o output.txt -k '3,2' -v '3' -n 'in1,in2,in3' input1.txt input2.txt input3.txt

END
        exit 1;
    }

    my @opt_input = @ARGV;

    my $key_specs       = parse_key_specs($opt_keys);
    my @key_columns     = map { $_->[0] } @$key_specs;
    my @val_columns     = parse_val($opt_columns); 
    my $comparator      = make_comparator($key_specs);
    my $files_and_nicks = parse_nicknames(\@opt_input, $opt_nicknames);

    # say Dumper \@opt_input, $files_and_nicks;

    for my $pair (@$files_and_nicks) { # foreach file and its nickname
        my ($input_file, $nick) = @$pair;

        # no option to disable like gff version-- can't trust user to sort properly. 
        sort_file_in_place($input_file, $key_specs); 

        # output to tempfile
        my (undef, $tempout_file) = tempfile("$opt_output.XXXXX");
        open my $tempout, '>', $tempout_file;

        my $input_parser = make_split_iterator($input_file);

        if (-f $opt_output){
            my $mashup_parser = make_split_iterator($opt_output);
            my $num_mashup_vals;
            while (defined(my $mashup_line = $mashup_parser->())){
                my @mashup_keys = @{$mashup_line}[0 .. scalar(@key_columns) -1]; 
                my @mashup_vals = @{$mashup_line}[scalar(@key_columns) .. $#$mashup_line];

                $num_mashup_vals //= @mashup_vals;

                # first (header) line?
                if ($mashup_parser->('numread') == 1){
                    say $tempout join "\t", @mashup_keys, @mashup_vals, map { "${nick}_$_" } @val_columns;
                }
                else{
                    if (! $input_parser->('eof')){ # input file done? 
                        INPUT:
                        while (defined(my $input_line = $input_parser->())){
                            my @input_keys = extract_columns($input_line, \@key_columns);
                            my @input_vals = extract_columns($input_line, \@val_columns);

                            my $cmp = $comparator->(\@input_keys, \@mashup_keys);

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

                                # put the input line back and get another mashup line. this is cleaner than
                                # trying to let mashup catchup here. 
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

        # no existing output file - create one. 
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

