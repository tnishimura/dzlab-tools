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

my $makensis = "/usr/bin/makensis";
my $perl     = "/usr/bin/perl";

my $bindir = $FindBin::Bin;
my $logfile;
my $outdir = $bindir;
my $pull;
#my $earliest = 1303945200; # april 27, 2011 16:00. don't build anything before this.

my $result = GetOptions (
    "logfile=s" => \$logfile,
    "outdir=s" => \$outdir,
    "pull"     => \$pull,
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

if ($pull){
    $logger->info("--pull");
    $logger->info(Git::Repository->run('pull'));
}

my $committime;
my $committime_str;
my $version;

# get committime
{
    my $iter = Git::Repository::Log::Iterator->new('-1');
    if (defined(my $log = $iter->next())){
        $committime     = $log->{author_gmtime}; 
        $committime_str = strftime("%Y.%m.%d-%H.%M",localtime($log->{author_gmtime})); 
        $logger->info("commit time = $committime_str");
    }
    else {
        $logger->logdie("can't get git-log?");
    }
}

# set version
{
    my ($desc) = Git::Repository->run(describe => '--tags');
    $logger->debug("raw describe output: $desc");
    if ($desc =~ /^\d+\.\d+\.\d+$/){ # HEAD is tagged with xx.xx.xx 
        $version = $desc;
        $logger->info("commit was tagged as $version");
    }
    else { # use timetamp
        $version = "CURRENT-$committime_str";
        $logger->info("commit not tagged, using timestamp $version");
    }
}

my $filename = "dzlab-tools-$version.exe";
my $dest     = catfile($outdir, $filename);

if (! -f $dest){
    #rm -rf help
    # pod2projdocs -o help -l . -except doc/index.pod -except '^tmp'
    $logger->info("creating $dest");
    system("$perl -wlpe 's/__VERSION__/$version/' installer.nsi.in > installer.nsi");
    system("$perl -wlpe 's/__VERSION__/$version/' dzlab-check.pl.in > dzlab-check.pl");
    $logger->info(`$makensis installer.nsi`);
    unlink "installer.nsi";
    unlink "dzlab-check.pl";
    move $filename, $dest;
}
else {
    $logger->info("$dest already exists, not remaking...");
}


=head1 NAME

_build_current.pl - build currently checked out dzlab-tools installer. If it's tagged, use tag name.
Other use timestamp.

=head1 SYNOPSIS

Usage examples:

 _build_current.pl --logfile some/where/log --outdir some/where/outdir;

=cut

