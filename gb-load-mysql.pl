#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use File::Which;

use FindBin;
use lib "$FindBin::Bin/lib";
use GBUtil;
use GFF::Statistics qw/gff_detect_width/;
use DZUtil qw/approximate_line_count/;

my ($user, $pass, $database, $host) = load_mysql_config();

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "create|c" => \my $create,
    "wiggle-dir|d=s" => \my $dir,
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result);  

if ($create){
    my $drop_stmt = "drop database $database; create database $database;";
    say STDERR $drop_stmt;
    system("mysql", "-u$user", "-p$pass", "-e", $drop_stmt);
}

my $bp_seqfeature_load = do{
    my $first = 1;
    my $exe = which("bp_seqfeature_load") // which("bp_seqfeature_load.pl");
    my @common_args = ("-u", $user,  "-p", $pass, "--dsn", "dbi:mysql:database=$database;host=$host", "--adaptor", 'DBI::mysql');

    sub {
        my $file = shift;
        my @cmd = ($exe, @common_args, ($create && $first ? ("--create") : ()), $file);
        say STDERR join " ", @cmd;
        system(@cmd);
        $first = 0;
    };
};


for my $file (@ARGV) {
    if ($file =~ /\.(fa|fasta|fas)$/i){
        say STDERR "preparing $file";
        my ($normalized, $meta) = prepare_fasta($file, $dir);
        $bp_seqfeature_load->($normalized);
        $bp_seqfeature_load->($meta);
    }
    elsif ($file =~ /\.gff3?$/i){
        if (approximate_line_count($file) > 1_000_000 && defined gff_detect_width($file)){
            say STDERR "converting $file to wig";
            gff_to_wig(
                file    => $file,
                dir     => $dir,
                compile => 1,
                bed     => 0,
                ctscore => 0,
            );
        }
        else{
            say STDERR "preparing $file";
            my $normalized = prepare_gff($file, $dir);
            $bp_seqfeature_load->($normalized);
        }
    }
}

