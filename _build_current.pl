#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use POSIX;
use File::Spec::Functions qw/catfile/;
use Git::Repository;
use Git::Repository::Log::Iterator;
use Pod::Usage;
use Getopt::Long;
use File::Copy;
use Config::General qw/ParseConfig/; 
use Pod::Usage;
use Getopt::Long;

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

chdir $FindBin::Bin;
my $r = Git::Repository->new(git_dir => $bindir . '/.git');

# $logger->info("pulling");
# my $pull = $r->run('pull');
my $iter = Git::Repository::Log::Iterator->new($r, '-1');

if (defined(my $log = $iter->next())){
    $logger->debug("most recent commit: " . Dumper $log);
    my $version  = "CURRENT-" . strftime("%Y.%m.%d-%H.%M",localtime($log->{author_gmtime}));
    my $filename = "dzlab-tools-" . "$version.exe";
    my $dest     = catfile($outdir, $filename);
    
    if (! -f $dest){
        #rm -rf help
        # pod2projdocs -o help -l . -except doc/index.pod -except '^tmp'
        $logger->info("creating $dest");
        system("perl -wlpe 's/__VERSION__/$version/' installer.nsi.in > installer.nsi");
        system("perl -wlpe 's/__VERSION__/$version/' dzlab-check.pl.in > dzlab-check.pl");
        system("makensis installer.nsi");
        unlink "installer.nsi";
        unlink "dzlab-check.pl";
        move $filename, $dest;
    }
    else {
        $logger->info("$dest already exists, not remaking...");
    }
}
