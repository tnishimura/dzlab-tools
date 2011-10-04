#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use List::Util qw/sum/;
use FindBin;
use lib "$FindBin::Bin/lib";
use lib "$FindBin::Bin/Tree-Range/lib";
use Ends::NeighborMapCollection;
use Scalar::Util qw/looks_like_number/;
use List::MoreUtils qw/all/;
use Pod::Usage;
use Getopt::Long;

my $result = GetOptions(
    'gff-annotation|g=s' => \my $gff_annotation,
    'bin-width|b=i'      => \(my $bin_width = 100),
    'distance|d=i'       => \(my $distance = 5000),
    'stop-flag|s=i'      => \(my $stop_flag = 2),
    'stop-distance|k=i'  => \(my $stop_distance = 1500),
    'three-prime|3'      => \(my $three_prime = 0),
    'five-prime|5'       => \(my $five_prime = 0),
    'extract-id|x=s'     => \(my $attribute_id = 'ID'),
);
if (! $result 
    || ! $gff_annotation
    || ! -f $gff_annotation 
    || ! ($five_prime xor $three_prime) 
    ){
    pod2usage(-verbose => 99);
}

my $nmc = Ends::NeighborMapCollection::new_cached(
    file            => $gff_annotation,
    tag             => $attribute_id,
    flag            => $stop_flag,
    distance        => $distance,
    prime           => $three_prime ? 3 : 5,
    flag_6_distance => $stop_distance,
    binwidth        => $bin_width,
);

=head1 USAGE

=head1 OPTIONS

=over

=item -g <annotation_gff> | --gff-annotation <annotation_gff>

=item -5 | --five-prime

=item -3 | --three-prime

=item  -b <width> | --bin-width <width>

=item  -d <distance> | --distance <distance>

=item  -x <id_tag> | --extract-id <id_tag>

=item -s | --stop-flag

=item  -k <distance> | --stop-distance <distance>

=back

=cut
