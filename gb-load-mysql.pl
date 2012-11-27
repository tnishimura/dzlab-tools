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

use GFF::Statistics qw/gff_detect_width/;
use DZUtil qw/approximate_line_count/;

my ($user, $pass, $database, $host) = load_mysql_config();

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "create|c" => \my $create,
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
    my @common_args = (
        "-u", $user,  
        "-p", $pass, 
        "--dsn", "dbi:mysql:database=$database;host=$host", 
        # "--adaptor", 'DBI::mysql',
    );

    sub {
        my $file = shift;
        my @cmd = ($exe, @common_args, ($create && $first ? ("--create") : ()), $file);
        say Dumper \@cmd;
        say STDERR join " ", @cmd;
        if (0 != system(@cmd)){
            die "failed: $?";
        }
        $first = 0;
        say "meow";
    };

};

my $config_file = shift // die "need config file";
my $config = LoadFile($config_file);

for my $fasta (@{$config->{fasta}}) {
    $bp_seqfeature_load->($fasta->{meta});
    # $bp_seqfeature_load->($fasta->{staging});
}

for my $gff (@{$config->{gff}}) {
    $bp_seqfeature_load->($gff->{staging});
}

for my $gffwig (@{$config->{gffwig}}) {
    $bp_seqfeature_load->($gffwig->{meta});
}



# --- sample 
# fasta:
#   - meta: /home/toshiro/demeter/staging/TAIR_reference.fas.normalized.meta
#     source: TAIR8
#     staging: /home/toshiro/demeter/staging/TAIR_reference.fas.normalized
# gff:
#   - feature:
#       - CDS
#       - chromosome
#       - exon
#       - five_prime_UTR
#       - gene
#       - mRNA
#       - mRNA_TE_gene
#     source: TAIR8
#     staging: /home/toshiro/demeter/staging/TAIR8_gmod.gff.normalized
# gffwig:
#  - feature: CHG
#    meta: /home/toshiro/demeter/staging/all.chg-col.w50-CHG-methyl.meta.gff
#    source: at-en-lerfie-x-col-wt-chg
#    type: methyl
#  - feature: CHG
#    meta: /home/toshiro/demeter/staging/all.chg-col.w50-CHG-coverage.meta.gff
#    source: at-en-lerfie-x-col-wt-chg
#    type: coverage
#  - feature: CG
#    meta: /home/toshiro/demeter/staging/all.cg-col.w50-CG-methyl.meta.gff
#    source: at-en-lerfie-x-col-wt-cg
#    type: methyl
#  - feature: CG
#    meta: /home/toshiro/demeter/staging/all.cg-col.w50-CG-coverage.meta.gff
#    source: at-en-lerfie-x-col-wt-cg
#    type: coverage
#  - feature: CHH
#    meta: /home/toshiro/demeter/staging/all.chh-col.w50-CHH-methyl.meta.gff
#    source: at-en-lerfie-x-col-wt-chh
#    type: methyl
#  - feature: CHH
#    meta: /home/toshiro/demeter/staging/all.chh-col.w50-CHH-coverage.meta.gff
#    source: at-en-lerfie-x-col-wt-chh
#    type: coverage
