#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Pod::Usage;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/lib";
use DZUtil qw/split_names/; 

my $output;
my $append;
my $help;
my $label;
my $inplace;
my $min = 7;
my $result = GetOptions (
    "append|a"    => \$append,
    "output|o=s"  => \$output,
    "label|l=s"   => \$label,
    "help|h"      => \$help,
    "minchar|c=i" => \$min,
    "inplace"     => \$inplace,
);

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $help || !@ARGV;

if ($output && @ARGV != 1){
    say "if --output (-o) is given, only 1 input may be given.";
    exit 1;
}

sub outname{
    my $input = shift;
    if ($inplace){
        return $input;
    }
    return split_names($input, 'renamed');
}

sub label{
    if ($label){
        return $label;
    }
    my $prefix = shift;
    $prefix =~ s{^[./]*}{};

    if ($prefix =~ m/(.*?)\.[^.]+$/){
        $prefix = $1;
    } 

    if (! $min){
        return $prefix;
    }
    elsif (length $prefix <= $min){
        return $prefix;
    } else{
        return substr $prefix, 0, $min;
    }
}

# create input => output map

my %files;
if ($output){ # single file
    $files{$ARGV[0]} = $output;
} 
else{
    %files = map {
        my @input = -f && $_ || glob("$_/*.gff");
        map {
            $_ => outname($_)
        } @input
    } @ARGV;
}

while (my ($in,$out) = each %files) {
    my $feature = label($in);
    my ($ifh, $ofh);

    if ($in eq $out){
        # inplace
        open $ifh, '<', $in;
        unlink $in;
        open $ofh, '>', $in;
    }
    else{
        open $ifh, '<', $in;
        open $ofh, '>', $out;
    }

    while (defined(my $line = <$ifh>)){
        chomp $line;
        my @s = split /\t/, $line;
        if ($append){
            $s[2] .= $feature;
        }
        else{
            $s[2] = $feature;
        }
        say $ofh join "\t", @s;
    }

    close $ifh;
    close $ofh;
}

=head1 NAME
 
feature_rename.pl - rename the feature column (col 3) in a GFF file, so it is differentiable in 
utilities like SignalMap.
 
=head1 SYNOPSIS

Read in input1.gff, rename the 3rd column to the first 7 letters (or less) of
the file name, create a file named output.gff

 feature_rename.pl -o output.gff input.gff

Like above, but output file names are assumed to be input1-renamed.gff, etc.
More than one input file allowed.

 feature_rename.pl input1.gff input2.gff ...

Rename everything in current directory ('.' means current directory).

 feature_rename.pl . 

Like above, but rename the 3rd column to 'myfeaturename' or whatever you give
to --label.

 feature_rename.pl --label myfeaturename input1.gff

=head1 OPTIONS

=over

=item -o <file> | --output <file>

=item -l <label> | --label <label>

=item -c <c> | --minchar <c>

=item -i | --in-place

=item -h | --help

=back

=cut

#=item --in-place
#Like the two examples above, BUT replace the input file with the renamed one
#(original overwritten). Careful.

# feature_rename.pl -i input1.gff
