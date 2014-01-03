#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Cwd qw/getcwd/;
use File::Basename qw/basename dirname/;
use File::Path qw/make_path remove_tree/;
use File::Spec::Functions qw/rel2abs canonpath catdir catfile updir/;
use File::Copy;
use FindBin;
use lib "$FindBin::Bin/lib";
use DZUtil qw/open_maybe_compressed/;

END {close STDOUT}
$| = 1;

my $fasta = shift || die "usage: $0 fasta";
my $basename = basename $fasta, qw/.fa .fasta .FA .FASTA/;

my $fh = open_maybe_compressed $fasta;

my $current_seq;
my $current_fh;
while (defined(my $line = <$fh>)){
    chomp $line;
    if ($line =~ /^>(\w+)/){
        if ($current_seq){
            close $current_fh;
        }
        $current_seq = $1;
        open $current_fh, '>', "$basename.$current_seq.fa";
    }
    if ($current_fh){
        say $current_fh $line;
    }
}
if ($current_fh){
    close $current_fh;
}

