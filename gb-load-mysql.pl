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

my $load = which("bp_seqfeature_load") // which("bp_seqfeature_load.pl");
my ($user, $pass, $database, $host) = load_mysql_config();

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "create|c" => \my $create,
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result);  

my $drop_stmt = "drop database $database; create database $database;";
if ($create){
    say STDERR $drop_stmt;
    system("mysql", "-u$user", "-p$pass", "-e", $drop_stmt);
}
my $first = 1;
for my $file (@ARGV) {
    if ($file =~ /\.(fa|fasta|fas)$/i){
    }
    elsif ($file =~ /\.(gff)$/i){
    }
    run($load, 
        "-u", $user,  
        "-p", $pass, 
        ($create && $first ? ("--create") : ()), 
        "--dsn", 
        "dbi:mysql:database=$database;host=$host", 
        "--adaptor", 'DBI::mysql',
        $file
    );
    $first = 0;
}

sub run{
    my @cmd = @_;
    say STDERR join " ", @cmd;
    system(@cmd);
}
