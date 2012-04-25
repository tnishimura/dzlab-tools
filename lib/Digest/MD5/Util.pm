package Digest::MD5::Util;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Digest::file qw(digest_file_hex);
use File::Basename;
use File::Spec::Functions qw/catfile abs2rel rel2abs/;

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(md5_of_file md5_confirm md5_create_checksum);

sub md5_of_file{
    my $file = shift;
    croak "md5_of_file needs file argument" if (! defined $file);
    croak "Does not exist: $file" if (! -f $file);
    return digest_file_hex($file, "MD5");
}

sub md5_confirm{
    my $checksum_file = shift;
    croak "md5_of_file needs file argument" if (! defined $checksum_file);
    croak "Does not exist: $checksum_file" if (! -f $checksum_file);

    my (undef, $basedir) = fileparse($checksum_file);
    my %results;

    open my $fh, '<', $checksum_file;
    while (defined(my $line = <$fh>)){
        chomp $line;
        if ($line =~ /^([a-f0-9]{32})\s\s(.*)$/){
            my $relpath = $2;
            my $fullpath = catfile($basedir, $relpath);
            my $expected = $1;
            my $got = md5_of_file ($fullpath);
            $results{$fullpath} = $got eq $expected;
        }
        else{
            croak "malformed $line";
        }
    }
    close $fh;

    if (0 == grep { $_ == 0 } values %results){
        return  [keys %results];
    }
    else{
        return;
    }
}

sub md5_create_checksum{
    my $checksum_file = rel2abs shift;
    my @files = map { rel2abs $_ } @_;

    my (undef,$path) = fileparse($checksum_file);
    open my $out, '>', $checksum_file;

    for my $f (@files) {
        say $out md5_of_file($f) . "  " . abs2rel($f, $path);
    }
    close $out;
}

1;

