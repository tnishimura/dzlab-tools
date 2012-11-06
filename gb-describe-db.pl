#!/usr/bin/env perl
# ripped mostly from gmod wiki
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Getopt::Long;
use Bio::DB::SeqFeature::Store;
use FindBin;
use lib "$FindBin::Bin/lib";
use GBUtil;

my ($user, $pass, $database, $host) = load_mysql_config();

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "user|u=s"     => \$user,
    "pass|p=s"     => \$pass,
    "database|d=s" => \$database,
    "host|h=s"     => \$host,
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || ! $user || !$pass || !$database || !$host);  
 
my $db = Bio::DB::SeqFeature::Store-> new( 
    -adaptor => 'DBI::mysql',
    -dsn     => "dbi:mysql:$database:$host",
    -user    => $user,
    -pass    => $pass,
    );

## Are sub-features being indexed?
print "\nIndexed sub-features? : ",
  $db->index_subfeatures, "\n";
 
## What serializer is being used?
print "\nThe serializer is : ",
  $db->serializer, "\n";
 
## List all the feature types in the database:
print "\nFeature types:\n";
print "\t", join("\n\t", $db->types), "\n";
 
## List how many there are for each feature type in the database
print "Feature type counts:\n";
my %types = $db->types(-count => 1);
print join("\n", map { $_ . "\t". "$types{$_}" } keys %types), "\n";
 
## List the feature attributes (tags) in the database:
print "\nAttributes:\n";
print "\t", join("\n\t", $db->attributes), "\n";
 
## How many sequence ids in the database:
print "\nThere are : ", scalar($db->seq_ids), " sequences in the database\n";
 
### How many features in the database:
print "There are : ", scalar($db->features), " features in the database\n";
