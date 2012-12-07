#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Storable qw/dclone/;
use YAML qw/Load Dump LoadFile DumpFile/;
use File::Which;

use FindBin;
use lib "$FindBin::Bin/lib";
use GBUtil;
use GBUtil;
use GBUtil::InputFile::MethylGFF;
use GBUtil::InputFile::GFF;
use GBUtil::InputFile::Fasta;

my ($user, $pass, $database, $host) = load_mysql_config();

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "create|c"     => \my $create,
    "user|u=s"     => \$user,
    "pass|p=s"     => \$pass,
    "database|d=s" => \$database,
    "host|h=s"     => \$host,
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result);  

if ($create){
    my $drop_stmt = "drop database $database; create database $database;";
    say STDERR $drop_stmt;
    system("mysql", "-u$user", "-p$pass", "-e", $drop_stmt);
}

sub bp_seqfeature_load{
    state $first = 1;
    state $exe = which("bp_seqfeature_load") // which("bp_seqfeature_load.pl");
    state $common_args = [
        "-u", $user,  
        "-p", $pass, 
        "--dsn", "dbi:mysql:database=$database;host=$host", 
        # "--adaptor", 'DBI::mysql',
    ];

    my $file = shift;
    my @cmd = ($exe, @$common_args, ($create && $first ? ("--create") : ()), $file);
    if (0 != system(@cmd)){
        die "failed: $?";
    }
    $first = 0;
}


my $config_file = shift // die "need config file";
my @input_files = @{LoadFile($config_file)};

for my $input (@input_files) {
    say STDERR $input;
    for my $f ($input->upload_files) {
        bp_seqfeature_load( $f );
    }
}
