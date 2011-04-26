#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Pod::Usage;
use Getopt::Long;
use Git::Repository;
use FindBin;
use File::Spec::Functions;

my $bindir = $FindBin::Bin;
my $logfile;
my $outdir = $bindir;

my $result = GetOptions (
    "logfile=s" => \$logfile,
    "outdir=s" => \$outdir,
);

if (!$result){
    say "usage: $0 --logfile some/where/log --outdir some/where/outdir";
    exit 1;
}

use Log::Log4perl qw/:easy/;
Log::Log4perl->easy_init( { 
    level    => $DEBUG,
    layout   => '%d{HH:mm:ss} %p> (%L) %M - %m%n',
    ($logfile ?  (file => ">>$logfile") : ()),
} );
my $logger = get_logger();

chdir $bindir;

$logger->info("pull: " . Git::Repository->run('pull'));
my @tags = Git::Repository->run('tag');

say join ",", @tags;

for my $tag (@tags) {
    if (! -f catfile($outdir,"dzlab-tools-$tag.exe")){
        $logger->info("$tag needs building");
        $logger->info("checkout: " . Git::Repository->run(checkout => $tag));
        system("./_build.pl $tag");
    } else{
        $logger->info("$tag already built, skipping");
    }
}

$logger->info("checkout master: " . Git::Repository->run(checkout => 'master'));
