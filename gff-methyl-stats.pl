#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Statistics qw/methylation_stats/;
use File::Spec::Functions qw/rel2abs/;
use YAML qw/LoadFile DumpFile/;
use Cwd qw/getcwd/;
use File::Basename qw/basename dirname/;
use File::Path qw/make_path remove_tree/;
use File::Spec::Functions qw/canonpath catdir catfile updir/;
use Parallel::ForkManager;

END {close STDOUT}
$| = 1;
use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "tmp-dir|d=s" => \my $tmpdir,
    "parallel|p=i" => \(my $parallel=1),
);
pod2usage(-verbose => 99) if (!$result || !@ARGV);  

my $pm = Parallel::ForkManager->new($parallel);

if (! $tmpdir || ! -d $tmpdir){
    $tmpdir = getcwd();
}


sub memoname{
    my $name = rel2abs shift;
    $name =~ s/\W/_/g;
    $name =~ s/^_//g;
    return catfile($tmpdir, $name);
}

my %all_stats;

for my $file (@ARGV) {
    $pm->start and next;
    my $memo = memoname($file);
    if (! -f $memo){
        my $stats = methylation_stats($file);
        DumpFile($memo, $stats);
        $all_stats{$file} = $stats;
    }
    $pm->finish; 
}
$pm->wait_all_children;

for my $file (@ARGV) {
    my $memo = memoname($file);
    if (-f $memo){
        $all_stats{$file} = LoadFile($memo);
    }
    else{
        die "why doesn't $memo exist?";
    }
}

# header

say join "\t", "", @GFF::Statistics::rownames;

for my $file (sort keys %all_stats) {
    say join "\t", $file, @{$all_stats{$file}}{@GFF::Statistics::rownames};
}


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
