#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;

use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) if ! $opt_database;

use DBI;

my $sth;
my $dbh = DBI->connect("dbi:SQLite:dbname=$opt_database","","", {RaiseError => 1, AutoCommit => 1});
$dbh->do("PRAGMA automatic_index = OFF");
$dbh->do("PRAGMA journal_mode = OFF");
$dbh->do("PRAGMA cache_size = 80000");

$dbh->do("create table scores (chr, coord numeric)");
$dbh->do("create table columns (filename, nickname, num)");
$dbh->do("create index idx1 on scores (chr,coord)");
$dbh->disconnect;

=head1 NAME
 
single_c_mashup_create.pl
 
=head1 SYNOPSIS

 single_c_mashup_create.pl

=head1 OPTIONS

=over

=item  -d <file> | --database <file>

=item -h | --help

=back

=cut

