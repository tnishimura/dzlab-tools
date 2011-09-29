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
use Storable;

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
    'zero'               => \(my $zero_flag_region = 0),
    'debug'              => \(my $debug = 0),
    'singleton|1'        => \(my $singleton = 0),
    #'no-skip'            => \(my $noskip),
    #'verbose|v'          => sub { use diagnostics; },
    #'quiet|q'            => sub { no warnings; },
    #'help|h'             => sub { pod2usage( -verbose => 1 ); },
    #'manual|m'           => sub { pod2usage( -verbose => 2 ); }
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
    \$output           = $output
    \$debug            = $debug
    \$zero_flag_region = $zero_flag_region
    \$singleton        = $singleton
LOGMSG

my $nmc = Ends::NeighborMapCollection::new_cached(
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

        RESLOOP:
        for my $result ($nmc->lookup($seq, $start, $end)) {
            my ($id, $bin, $overlap) = ($result->{item}[0], $result->{item}[1], $result->{overlap});
            add_to_table($id, $bin, $c, $t, $overlap, $len);
            last RESLOOP if $singleton;
        }
    }
}

dump_table($output);
dump_average($output eq '-' ? '-' : $output . ".avg");

#######################################################################
{
    my %table; # { id => [[current avg,num], [1,2], ... ] }
    my $numbins;
    my $distance;
    my $binwidth;
    my @all_id;
    sub build_table{
        my ($nmc) = @_;
        $numbins = $nmc->numbins();
        $distance = $nmc->distance();
        $binwidth = $nmc->binwidth();
        @all_id = $nmc->all_id();

        for my $id ($nmc->all_id) {
            $table{$id} = [map {$zero_flag_region && $nmc->bin_valid($id, $_) ? [0,0] : undef} (0 .. $numbins - 1)];
        }
    }
    sub add_to_table{
        my ($id, $bin, $c, $t, $overlap, $len) = @_;
        #die " $id, $bin, $c, $t, $overlap, $len";
        my $score = safediv($c, $c+$t) * $overlap / $len;
        if (! defined $table{$id}[$bin]){
            $table{$id}[$bin] = [$score,1];
        }
        my ($current, $num) = @{$table{$id}[$bin]};
        $table{$id}[$bin][0] = ($current * $num + $score) / ($num + 1);
        $table{$id}[$bin][1] = $num + 1;
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

        for my $id (sort @all_id) {
            say $fh join "\t", 
            $id, 
            map { 
                defined $_ ? $_->[0] : 'na'
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
            my @scores = 
            map {$_->[0]} 
            grep {defined} 
            map { $table{$_}[$bin] } 
            @id_list;

            say $fh (-$distance + $bin * $binwidth), "\t", safediv(sum(@scores), scalar(@scores));
        }
    }
}

=head1 NAME

armageddon_analysis.pl - Given a GFF annotation, create a histogram of scores
relative to ends of annotation entries from input GFF files.  Can also be called as arma.pl

=head1 USAGE

Note: arma.pl is an alias for armageddon_analysis.pl

At a bare minimum, you need to specify an annotation (-g), an end (-5 or -3), an input, and output (-o). Remember that stop flag (-s) is 2 by default.

 armageddon_analysis.pl -g annotation.gff -o output.ends -5 input.gff
 armageddon_analysis.pl -g annotation.gff -o output.ends -3 input.gff
 arma.pl -g annotation.gff -o output.ends -5 input.gff
 arma.pl -g annotation.gff -o output.ends -3 input.gff

Stop flag 6 analysis with stop-distance (-k) of 150:

 arma.pl -g annotation.gff -o output.ends -5 -s 6 -k 150 input.gff

You can fiddle with bin width (-b) and absolute distance (-d) too, though I don't think anybody does this.

 arma.pl -g annotation.gff -b 100 -d 5000 -s 2 -5 -o output.ends input.gff

For non-arabidopsis annotations, you'll have to change the extract ID (-x). For rice, it's Alias (I think):

 arma.pl -g annotation.gff -5 -x Alias input.gff

Emulate (technically incorrect) ends_analysis.pl behavior with --singleton (see Differences section below):

 arma.pl --singleton -g anno.gff -5 -s 2 -o output.ends input.gff

Zero all valid regions (--zero) without any mapping elements:

 arma.pl --zero -g anno.gff -5 -s 0 -o output.ends input.gff
 arma.pl --zero -g anno.gff -5 -s 2 -o output.ends input.gff

=head1 OPTIONS

Most of the options are the same as for ends_analysis.pl.

=over

=item -g <annotation_gff> | --gff-annotation <annotation_gff>

=item -o <file> | --output <file>

Output file name.  Average ends file will be the same name suffixed with ".avg".  (something like "outputfile.ends.avg").

=item -5 | --five-prime

Center analysis on 5' end.

=item -3 | --three-prime

Center analysis on 3' end.

=item  -b <width> | --bin-width <width>

Histogram bin width.  Default 100.  

=for Euclid
    width.default:     100

=item  -d <distance> | --distance <distance>

distance (in base pairs) from end terminal to search, both ways. Default 5000.

=for Euclid
    distance.default:     5000


=item  -x <id_tag> | --extract-id <id_tag>

Annotation's column 9 ID name.  Each annotation line needs a unique name, and it will use this option's value to extract from column 9.  For example, for TAIR, column 9 will have something like "ID=AT1G12345", in which case you should use the (default) value of 'ID'.

=for Euclid
    id_tag.default:     'ID'

=item  --zero 

If a bin has no scores mapping to it, but it is valid for analysis according to stop flag, output 0 instead of 'na'.

=item -s | --stop-flag

Defines the range of a particular annotation entry's analysis.  Can be 0, 2, or 6.  See Stop Flag section below. Default 2.

=item  -k <distance> | --stop-distance <distance>

Distance from genes to stop from (stop flag 6 only, see Stop Flag section.) Default 15000.

=for Euclid
    distance.default:     1500

=back

=head2 Flag options

=over

=item -s 0

Don't look at adjacent genes, but just keep going the maximun
distance away.  So, when  

=item -s 2

Look within the gene, but look upstream or downstream in intergenic
regions (stop at adjacent gene).  So, when analysis is centered 5', will stop
at the 3' end, or at an adjacent gene, whichever comes first, but will not
include intergenic regions downstream.  It will include intergenic regions
upstream of the 5' end.  

=item -s 6

The same as stop flag #2, but stopping alignment by the
user-specified amount before the end of the gene.

=back

=head2 Differences vs. ends_analysis.pl

In the list below, the original ends_analysis.pl is refered to as Ends, and armageddon_analysis.pl as Armageddon.

=over 

=item armageddon_analysis.pl can handle uneven window sizes

Armageddon, unlike Ends, will correctly handle uneven window sizes. This means that you can pass arbitrary GFF files as input, not just single-c's.  In the situation where a single input GFF entry overlaps with multiple bins, its score is split according to overlap.  For example, if you have an input GFF entry with score .9 from 170 to 230 (length 61), and it overlaps with bins from 100-199 and 200-230, the first bin would get a contribution of .9 * 30 / 61 and the second would get .9 * 31/61.

=item armageddon_analysis.pl is more acurate

Several bugs were discovered in Ends during the development of Armageddon, and were rectified.  In particular, in Ends, if a single-c input position mapped to two bins of different genes (or transposons or whatever), it randomly chose one. 

For example, on a 5' analysis, if there were a gene from 5000 to 6000 on the 3'->5' strand, and another gene from 7000 to 8000 on the 5'->3', any position between 6001 and 6999 should influence both genes' analyses.  However, Ends randomly chose only one gene.  Armageddon handles this situation correctly.  As a result, Armageddon bin averages tend to be a tad higher than Ends.  For backwards compatibility with Ends, however, you can pass the --singleton flag to emulate Ends behavior in Armageddon.

In addition, Armageddon correctly ignores and genes that have their 5' or 3' end completely overlapped by another gene, unlike Ends.  Armageddon also fixes other severe bug rare bugs.

=item armageddon_analysis.pl automatically creates an average ends files for you

Armageddon will automatically create an averaged ends file with the *.avg suffix appended to the output file (so if you passed -o output.ends, an average file called output.ends.avg will be created automatically).  This file only contains averages only (unlike average_ends_new.pl output which produces standard deviation, quartiles, and other columns no one seems to use).

=item armageddon_analysis.pl will create an annotation 'cache-file'

(This is actually a disadvantage). The first time you run Armageddon for a particular set of options, it will create a 'cache-file' to index and pre-process the annotation file.  These will be created in the same directory as the annotation file.  This is a slow process, so the first time, expect it to spend 3-10 minutes creating this file.  Subsequent runs will reuse the cache file, so it should be about as fast as Ends, though there is room for improvement. Note that it will create separate cache files for different options, ie, it will create two files if you run a -s 0 and a -s 2 analysis.  

=item armageddon_analysis.pl will compute scores from the 'c' and 't' values from column 9

Ends will trust column 6, where Armageddon will calculate the methylation score from the c/(c+t).  I may add an option to use column 6 directly if people have a problem.

=back 

=cut

#stop==1: look at adjacent genes and stop looking only when you encounter
#another gene. ie. does NOT stop at the 3' end when aligning at the 5' end,
#but keeps going into intergenic regions until it hits another gene. 
#
#stop==3: only look within a gene, so stop at the end of the gene.
#
#stop==4: start at the start of the gene if there's an overlap, otherwise go
#upstream.  stop at the end of the gene. 
#
#stop==5: start the given distance upstream of the gene end (ignore adjacent
#or overlapping genes), and stop at the end of the gene.  So if aligning at
#the 5' end, will start x bp upstream and stop at the 3' end.
