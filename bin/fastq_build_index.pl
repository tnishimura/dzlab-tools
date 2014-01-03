#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use DBI;
use Log::Log4perl qw/:easy/;
use File::Temp qw/mktemp/;

Log::Log4perl->easy_init( { 
    level    => $DEBUG,
    layout   => '%d{HH:mm:ss} %p> (%L) %M - %m%n',
} );
my $logger = get_logger();

die "usage: $0 sequence.fastq" unless @ARGV==1;
open my $fh, '<', $ARGV[0];

my $outfile = $ARGV[0] . ".index";
unlink $outfile if -f $outfile;

my $tmpfile = mktemp($outfile . ".tmp.XXXXX");

#######################################################################
# Get file size

my $size = (stat($ARGV[0]))[7];

#######################################################################
# Config db

my $sth;
my $dbh = DBI->connect("dbi:SQLite:dbname=$tmpfile","","", {RaiseError => 1, AutoCommit => 1});
$dbh->do("PRAGMA automatic_index = OFF");
$dbh->do("PRAGMA journal_mode = OFF");
$dbh->do("PRAGMA cache_size = 80000");
$dbh->do("create table fastq (id,sequence)");
$sth = $dbh->prepare("insert into fastq (id,sequence) values (?,?)");
$dbh->{AutoCommit} = 0;

#######################################################################
# Loop

my $progress_step       = .05;
my $progress_checkpoint = .05;
my $commit_size         = 20000;
my $counter             = 0;
my $position            = tell $fh;

while (defined(my $id = <$fh>)){
    chomp $id;
    $logger->logdie("malformatted line $id") unless $id =~ /^@/;

    $id =~ s/^@//;


    my $sequence;
    if (! defined($sequence = <$fh>)){
        $logger->logdie("number of lines not a multiple of 4"); 
    }
    chomp $sequence;

    $sth->execute($id,$sequence);
    if (++$counter % $commit_size == 0){
        #$logger->debug("$counter");
        $dbh->commit or $logger->logdie("Couldn't commit?");
    }
    if ($position / $size > $progress_checkpoint){
        $logger->debug(sprintf("%.2f%%" , 100 * $position / $size));
        $progress_checkpoint += $progress_step;
    }

    # skip quality lines
    for (1..2){
        if (! defined <$fh>){ $logger->logdie("number of lines not a multiple of 4"); }
    }

    $position = tell $fh;
}
$logger->debug("Done inserting. Creating index");

$dbh->commit;
$dbh->{AutoCommit} = 1;
$dbh->do("create index idx1 on fastq (id)");

$logger->debug("Done creating index");

$dbh->disconnect or $logger->logdie("Couldn't disconnect?");

rename $tmpfile, $outfile or $logger->logdie("Couldn't create $outfile");
