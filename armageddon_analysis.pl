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
use List::MoreUtils qw/all/;
use Pod::Usage;
use Getopt::Long;
use GFF::Parser;
use DZUtil qw/safediv/;
use Counter;

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
#$gff_annotation = "/wip/tools/annotations/AT/gmod/TAIR8_genes.gff";
#$five_prime = 1;
if (! $result 
    || ! $gff_annotation
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
build_table($nmc);
say STDERR "Table built";
my $counter = Counter->new(verbose => 1);

for my $file (@ARGV) {
    my $parser = GFF::Parser->new(file => $file);
    while (defined(my $gff = $parser->next())){
        $counter->increment();
        my ($seq, $start, $end, $c, $t) = 
        ($gff->sequence(), $gff->start(), $gff->end(), $gff->get_column('c'), $gff->get_column('t'));
        next unless all { defined($_) } ($seq, $start, $end, $c, $t);
        my $len = $end-$start+1;

        for my $result ($nmc->lookup($seq, $start, $end)) {
            my ($id, $bin, $overlap) = ($result->{item}[0], $result->{item}[1], $result->{overlap});
            add_to_table($id, $bin, $c * $overlap / $len, $t * $overlap / $len);
        }
    }
}

dump_table($output);
dump_average($output eq '-' ? '-' : $output . ".avg");

#######################################################################
{
    my %table; # { id => [[0,0], [1,2], ... ] }
    my $numbins;
    my $distance;
    my $binwidth;
    sub build_table{
        my ($nmc) = @_;
        $numbins = $nmc->numbins();
        $distance = $nmc->distance();
        $binwidth = $nmc->binwidth();

        %table = map { 
            $_ => [map {undef} (1 .. $numbins )] 
        }
        $nmc->all_id;
    }
    sub add_to_table{
        my ($id, $bin, $c, $t) = @_;
        if (! defined $table{$id}[$bin]){
            $table{$id}[$bin] = [0,0];
        }
        $table{$id}[$bin][0]+=$c;
        $table{$id}[$bin][1]+=$t;
    }
    sub dump_table{
        my $file = shift;
        my $fh;
        if ($file eq '-'){
            $fh = \*STDOUT;
        }
        else {
            open $fh, '>', $file;
        }

        for my $id (sort keys %table) {
            say $fh join "\t", 
            $id, 
            map { 
                if (defined){
                    my ($c, $t) = @$_;
                    $c + $t == 0 ? 'na' : $c/($c+$t);
                }
                else{
                    'na';
                }
            } @{$table{$id}};
        }
    }
    sub dump_average{
        my $file = shift;
        my $fh;
        if ($file eq '-'){
            $fh = \*STDOUT;
        }
        else {
            open $fh, '>', $file;
        }

        my @id_list = keys %table;
        for my $bin (0 .. $numbins - 1) {
            my @scores = map {safediv($_->[0], $_->[0] + $_->[1])} grep {defined} map { $table{$_}[$bin] } @id_list;

            say $fh (-$distance + $bin * $binwidth), "\t", safediv(sum(@scores), scalar(@scores));
        }
    }
}
