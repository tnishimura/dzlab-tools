#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Scalar::Util qw/looks_like_number/;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Sam::Parser; 
use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    #"fix-rc|rc"       => \(my $fixrc),
    "output|o=s" => \(my $output = '-'),
);
if (! $result || (! @ARGV && -t STDIN)){
    pod2usage(-verbose => 2, -noperldoc => 1);
}
my $outfh = $output eq '-' ? \* STDOUT : IO::File->new($output, 'w');

#######################################################################
# now produce stats from intermediate files

use Hash::Util qw/lock_keys unlock_keys lock_hash unlock_hash/;

my %stats = (
    num_reads                   => 0,
    num_aligned_concordantly    => 0,
    num_aligned_disconcordantly => 0,
);
lock_keys(%stats);

my $previous_id;
my $lastid_already_counted;
my $total_bases = 0;

my $parser = Sam::Parser->new(file => \*ARGV, skip_unmapped => 0);
while (defined(my $sam1 = $parser->next)){
    # say "$.: " .  $sam1->readid;
    # say "- $. read: " . $sam1->readid;
    my $sam2 = $parser->next;
    if (! $sam2){
        die "uneven number of reads?";
    }
    # printf("%s and\n%s\n", $sam1->readid, $sam2->readid);
    if ($sam1->readid ne $sam2->readid){
        die sprintf("reads not in pairs?\n%s and\n%s", $sam1->readid, $sam2->readid);
    }
    my $concordant = (
        $sam1->is_mapped && $sam2->is_mapped &&
        (
            ( $sam1->is_first_segment && $sam2->is_second_segment ) xor
            ( $sam2->is_first_segment && $sam1->is_second_segment )
        )
    );
    my $disconcordant = (
        (! $sam1->is_mapped &&   $sam2->is_mapped ) xor 
        (  $sam1->is_mapped && ! $sam2->is_mapped ) 
    );
    $stats{num_reads}++;
    $stats{num_aligned_concordantly}++ if $concordant;
    $stats{num_aligned_disconcordantly}++ if $disconcordant;

    say STDERR $stats{num_reads} if $stats{num_reads} % 10000 == 0;

}

my $concordant_percentage    = $stats{num_reads} > 0 ? sprintf("%0.2f", 100 * $stats{num_aligned_concordantly} / $stats{num_reads}) : "0";
my $disconcordant_percentage = $stats{num_reads} > 0 ? sprintf("%0.2f", 100 * $stats{num_aligned_disconcordantly} / $stats{num_reads}) : "0";

$outfh->print(<<END);
number of reads total: $stats{num_reads}
number of reads aligning concordantly: $stats{num_aligned_concordantly} ($concordant_percentage%)
number of reads aligning disconcordantly: $stats{num_aligned_disconcordantly} ($disconcordant_percentage%)
END

$outfh->close;

=head1 NAME

Assumes that there is only 1 alignment per (which is the default mode in
bowtie2).

 sam-sam-pe.pl bowtie2-paired-ends-alignment.sam 

=cut


