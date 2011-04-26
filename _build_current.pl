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

use Log::Log4perl qw/:easy/;
Log::Log4perl->easy_init( { 
    level    => $DEBUG,
    #file     => ">run.log",
    layout   => '%d{HH:mm:ss} %p> (%L) %M - %m%n',
} );
my $logger = get_logger();

my $output = $FindBin::Bin;
my $result = GetOptions ( "output|o=s" => \$output,);
if (!$result || !$output){
    say "usage: $0 -o output_dir";
    exit 128;
}

chdir $FindBin::Bin;
my $r = Git::Repository->new(git_dir => $FindBin::Bin . '/.git');

# $logger->info("pulling");
# my $pull = $r->run('pull');
my $iter = Git::Repository::Log::Iterator->new($r, '-1');

if (defined(my $log = $iter->next())){
    $logger->debug("most recent commit: " . Dumper $log);
    my $version  = "CURRENT-" . strftime("%Y.%m.%d-%H.%M",localtime($log->{author_gmtime}));
    my $filename = "dzlab-tools-" . "$version.exe";
    my $dest     = catfile($output, $filename);
    
    if (! -f $dest){
        #rm -rf help
        # pod2projdocs -o help -l . -except doc/index.pod -except '^tmp'
        $logger->info("creating $dest");
        $logger->info(`perl -wlpe 's/__VERSION__/$version/' installer.nsi.in > installer.nsi"`);
        $logger->info(`perl -wlpe 's/__VERSION__/$version/' dzlab-check.pl.in > dzlab-check.pl`);
        $logger->info(`makensis installer.nsi`);
        move $filename, $dest;
    }
    else {
        $logger->info("$dest already exists, not remaking...");
    }
}


__DATA__
archive = /var/www/lighttpd/dzlab-tools/archive
www = /var/www/lighttpd/dzlab-tools
repo = /home/vaquero/dzlab-tools/.git
