#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use DBI;
use Log::Log4perl qw/:easy/;

Log::Log4perl->easy_init( { 
    level    => $DEBUG,
    layout   => '%d{HH:mm:ss} %p> (%L) %M - %m%n',
} );
my $logger = get_logger();


die "usage: $0 sequence.fastq" unless @ARGV==1;
open my $fh, '<', $ARGV[0];

my $outfile = $ARGV[0] . ".index";
unlink $outfile if -f $outfile;

my $sth;
my $dbh = DBI->connect("dbi:SQLite:dbname=$outfile","","", {RaiseError => 1, AutoCommit => 1});
$dbh->do("PRAGMA automatic_index = OFF");
$dbh->do("PRAGMA journal_mode = OFF");
$dbh->do("PRAGMA cache_size = 80000");
$dbh->do("create table fastq (id,byte)");
$sth = $dbh->prepare("insert into fastq (id,byte) values (?,?)");
$dbh->{AutoCommit} = 0;

my $counter = 0;
my $position = tell $fh;
while (defined(my $line = <$fh>)){
    chomp $line;
    die "malformatted line $line" unless $line =~ /^@/;
    $sth->execute($line,$position);
    if (++$counter % 20000 == 0){
        $logger->debug("$counter");
    }
    <$fh>;
    <$fh>;
    <$fh>;
    $position = tell $fh;
}
$dbh->commit;
$dbh->{AutoCommit} = 1;
$dbh->do("create index idx1 on fastq (id)");


