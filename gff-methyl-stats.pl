#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;

use Getopt::Long;
use Parallel::ForkManager;
use Pod::Usage;
use YAML qw/Dump/;

use FindBin;
use lib "$FindBin::Bin/lib";
use DZUtil qw/memodo memofile gimmetmpdir/;
use GFF::Statistics qw/methylation_stats/;

END {close STDOUT}
$| = 1;

my $result = GetOptions (
    "tmp-dir|d=s"  => \my $tmpdir,
    "parallel|p=i" => \(my $parallel=1),
    "force|f"      => \(my $force),
    "wp=s" => \(my $wanted_percentiles = ".05,.25,.50,.75,.95"),
);

pod2usage(-verbose => 99) if (!$result || !@ARGV);  

$tmpdir = gimmetmpdir($tmpdir);
my %all_stats;

my $pm = Parallel::ForkManager->new($parallel);
$pm->run_on_finish(sub{ 
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $ref) = @_;
        if (defined($ref)) {  
            my ($file, $stats) = %$ref;
            $all_stats{$file} = $stats;
        } 
    });

for my $file (@ARGV) {
    $pm->start and next;
    my $memo = memofile($file, $tmpdir);

    my $stats = memodo($memo, sub{
            my ($stats) = methylation_stats($file, split /,/, $wanted_percentiles);
            return $stats;
        }, $force);

    $pm->finish(0, {$file => $stats}); 
}
$pm->wait_all_children;

say Dump \%all_stats;

=head1 NAME

gff-stats.pl - produce statistics about gff files

=head1 SYNOPSIS

Usage examples:

 gff-stats.pl [-d /some/temporary/directory/] file1.gff file2.gff ...

will produce a row for each file with the following columns for each:

=head1 Output Columns

=over

=item nuc_ct_mean

=item chr_ct_mean

=item mit_ct_mean

Mean C+T score per methylation site.  (Calculate C+T for each methylation site
in the gff file, calculate the mean.)

=item nuc_ct_median

=item chr_ct_median

=item mit_ct_median

Median C+T score per methylation site.  (Calculate C+T for each methylation
site in the gff file, calculate the median.)

=item nuc_methyl_mean

=item chr_methyl_mean

=item mit_methyl_mean

Mean C/(C+T) score per methylation site.  (Calcualte C/(C+T) for each
methylation site in the gff file, calculate the mean.)

=item nuc_methyl_total

=item chr_methyl_total

=item mit_methyl_total

(Total C count in GFF file)/(Total C+T count in GFF file).  This is different
from the methyl_mean scores b/c this is in total.  methyl_mean gives the same
weight to a site where c=1;t=2 as c=10;t=20, whereas methyl_total combines the
c and t's first.  

=item coverage

=back

=head1 OPTIONS

=over

=item  -d <dir> | --tmp-dir <dir>

Temporary directory for cache.

=back

=cut
