package GFF::Split;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;
use GFF::Parser;
use File::Basename;
use File::Spec::Functions;
use FindBin;
use File::Temp qw/mktemp/;
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
    my $file = shift;
    return _gff_split($file,'feature');
}
sub split_sequence{
    my $file = shift;
    return _gff_split($file,'sequence');
}

sub _gff_split{
    my ($file,$col) = (@_);
    my %fh;
    my @files;
    my %tempfiles; # temp => real 

    my $p = GFF::Parser->new(file => $file);
    GFF:
    while (defined(my $gff = $p->next())){
        my $part;
        if ($col eq 'feature'){
            $part = $gff->feature;
        } elsif ($col eq 'sequence'){
            $part = $gff->sequence;
        } else{
            croak "\$col needs to be sequence or feature";
        }

        next GFF if ! defined $part;

        if (! exists $fh{$part}){
            my $split_file = _split_name($file, $part);
            my $temp  = mktemp($split_file . ".tmpXXXX");
            open $fh{$part}, '>', $temp;
            push @files, $temp;
            $tempfiles{$temp} = $split_file;
        }

        say {$fh{$part}} $gff->to_string;
    }

    for my $handle (values %fh) {
        close $handle;
    }
    while (my ($temp,$real) = each %tempfiles) {
        rename $temp, $real;
    }
    return sort values %tempfiles;
}

1;

