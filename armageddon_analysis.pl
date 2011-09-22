#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use List::Util qw/sum/;
use FindBin;
use lib "$FindBin::Bin/lib";
use Ends::NeighborMapCollection;
use Scalar::Util qw/looks_like_number/;
use Pod::Usage;
use Getopt::Long;

END {close STDOUT}
$| = 1;

my $result = GetOptions(
    'gff-annotation|g=s' => \my $gff_annotation,
    'bin-width|b=i'      => \(my $bin_width = 100),
    'distance|d=i'       => \(my $distance = 5000),
    'stop-flag|s=i'      => \(my $stop_flag = 2),
    'stop-distance|k=i'  => \(my $stop_distance = 1500),
    'three-prime|3'      => \(my $three_prime = 0),
    'five-prime|5'       => \(my $five_prime = 0),
    'extract-id|x=s'     => \(my $attribute_id = 'ID'),
    'output|o=s'         => \(my $output = '-'),
    'debug'              => \(my $debug = 0),
    #'zero'               => \$zero_flag_region,
    #'no-skip'            => \(my $noskip),
    #'verbose|v'          => sub { use diagnostics; },
    #'quiet|q'            => sub { no warnings; },
    #'help|h'             => sub { pod2usage( -verbose => 1 ); },
    #'manual|m'           => sub { pod2usage( -verbose => 2 ); }
);
$gff_annotation = "t/test1.gff";
$five_prime = 1;
if (! $result 
    || ! -f $gff_annotation 
    || ! ($five_prime xor $three_prime) 
    ){
    say "usage: armageddon_analysis.pl ...";
    exit 1;
}

print STDERR <<"LOGMSG";
    \$gff_annotation  = $gff_annotation  
    \$bin_width       = $bin_width
    \$distance        = $distance
    \$stop_flag       = $stop_flag
    \$stop_distance   = $stop_distance
    \$three_prime     = $three_prime
    \$five_prime      = $five_prime
    \$attribute_id    = $attribute_id
    \$output          = $output
    \$debug           = $debug
LOGMSG

my $nmc = Ends::NeighborMapCollection->new(
    file            => $gff_annotation,
    tag             => $attribute_id,
    flag            => $stop_flag,
    distance        => $distance,
    prime           => $three_prime ? 3 : 5,
    flag_6_distance => $stop_distance,
    binwidth        => $bin_width,
);
#build_table($nmc);
#add_to_table('AT1G01910', 10, 100,100);
#add_to_table('AT1G01910', 10, 110,200);
#dump_table();

#######################################################################
{
    my %table; # { id => [[0,0], [1,2], ... ] }
    sub build_table{
        my ($nmc) = @_;
        %table = map { 
            $_ => [map {[0,0]} (1 .. $nmc->numbins )] 
        }
        $nmc->all_id;
    }
    sub add_to_table{
        my ($id, $bin, $c, $t) = @_;
        $table{$id}[$bin][0]+=$c;
        $table{$id}[$bin][1]+=$t;
    }
    sub dump_table{
        for my $id (sort keys %table) {
            say join "\t", 
            $id, 
            map { 
                my ($c, $t) = @$_;
                $c + $t == 0 ? 'na' : $c/($c+$t);
            } @{$table{$id}};
        }
    }
}
