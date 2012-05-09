package GFF::Split;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;
use File::Basename;
use File::Spec::Functions;
use FindBin;
use File::Temp qw/mktemp/;
use List::MoreUtils qw/all/;
use lib "$FindBin::Bin/lib";

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(split_feature split_sequence split_names);
our @EXPORT = qw();

=head2 split_names("/home/test/path.txt", qw/chr1 chr2 chr3/)

Return files that would be produced by splitting.

=cut

sub split_names{
    my ($fullpath, @splits) = @_;
    return sort map { _split_name($fullpath, $_) } @splits;
}

sub _split_name{
    my ($fullpath, $part) = @_;
    my ($filename, $path) = fileparse($fullpath);

    my $ext = ($filename =~ s/\.([^.]+)$//) ? $1 : "";

    return catfile($path, $filename) . "-$part" . ($ext ? ".$ext" : "");
}

sub split_feature{
    my ($file, @groups) = @_;
    return _gff_split($file,'feature',@groups);
}
sub split_sequence{
    my ($file, @groups) = @_;
    return _gff_split($file,'sequence',@groups);
}

sub _gff_split{
    my ($file,$col,@groups) = (@_);
    my %part2fh;
    my %temp2real; # temp => real 
    my %part2real;

    if (@groups){
        my %expected = map { $_ => _split_name($file, $_) } @groups;
        if (all { -f } values %expected){
            return %expected;
        }
    }

    open my $input, '<', $file;
    GFF:
    while (defined(my $line = <$input>)){
        chomp $line;
        next GFF if $line =~ /^\s*#/;
        my @split = (split /\t/, $line);
        next GFF if ! @split == 9;
        my ($sequence, $feature) = @split[0,2];

        my $part = 
        $col eq 'feature' ? $feature : 
        $col eq 'sequence' ? $sequence :
        croak "\$col needs to be sequence or feature";

        next if ($part eq '.' || $part eq '');

        if (! exists $part2fh{$part}){
            my $split_file = _split_name($file, $part);
            my $temp  = mktemp($split_file . ".tmp.XXXXX");
            open $part2fh{$part}, '>', $temp;
            $temp2real{$temp} = $split_file;
            $part2real{$part} = $split_file;
        }

        say {$part2fh{$part}} $line;
    }

    for my $handle (values %part2fh) {
        close $handle;
    }
    while (my ($temp,$real) = each %temp2real) {
        rename $temp, $real;
    }
    return %part2real;
}

1;

