#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use List::Util qw/sum/;
use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../Tree-Range/lib";
use Ends::NeighborMapCollection;
use Scalar::Util qw/looks_like_number/;
use List::MoreUtils qw/all/;
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
    'flag-7-num-bins|n=i' => \(my $flag_7_numbins),
);
if (! $result 
    || ! $gff_annotation
    || ! -f $gff_annotation 
    || ! ($five_prime xor $three_prime) 
    ){
    pod2usage(-verbose => 99);
}

print STDERR <<"LOGMSG";
    \$gff_annotation   = $gff_annotation  
    \$bin_width        = $bin_width
    \$distance         = $distance
    \$stop_flag        = $stop_flag
    \$stop_distance    = $stop_distance
    \$three_prime      = $three_prime
    \$five_prime       = $five_prime
    \$attribute_id     = $attribute_id
LOGMSG
print STDERR "    \$flag_7_numbins   = $flag_7_numbins" if $flag_7_numbins;

my $nmc = Ends::NeighborMapCollection::new_cached(
    file            => $gff_annotation,
    tag             => $attribute_id,
    flag            => $stop_flag,
    distance        => $distance,
    prime           => $three_prime ? 3 : 5,
    flag_6_distance => $stop_distance,
    binwidth        => $bin_width,
    numbins         => $flag_7_numbins,
);

=head1 NAME

 arma_build.pl -g anno.gff [-5 | -3] -b 100 -d 5000 -x ID -s [0 | 2 | 6 | 7] 

 Only Flag 6: -k 1500
 Only Flag 7: -n 100

=cut
