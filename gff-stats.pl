#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Statistics;
use File::Spec::Functions qw/rel2abs/;
use YAML qw/LoadFile DumpFile/;
use Cwd qw/getcwd/;
use File::Basename qw/basename dirname/;
use File::Path qw/make_path remove_tree/;
use File::Spec::Functions qw/canonpath catdir catfile updir/;


END {close STDOUT}
$| = 1;


sub memoname{
    my $name = rel2abs shift;
    $name =~ s/\W/_/g;
    $name =~ s/^_//g;
    return $name;
}

#my $test_cg = "t/big-cg.gff";
#my $test_chg = "t/big-chg.gff";
#my @opt_files = ($test_cg, $test_chg);

my %all_stats;

for my $file (@ARGV) {
    my $memo = memoname($file);
    if (-f $memo){
        $all_stats{$file} = LoadFile($memo);
    }
    else{
        my $stats = getstats($file);
        DumpFile($memo, $stats);
        $all_stats{$file} = $stats;
    }
}

# header

my @rownames = qw/
coverage
nuc_ct_mean nuc_ct_median 
chr_ct_mean chr_ct_median 
mit_ct_mean mit_ct_median 
nuc_methyl_mean chr_methyl_mean mit_methyl_mean 
/;

say join "\t", "", @rownames;

for my $file (sort keys %all_stats) {
    say join "\t", $file, @{$all_stats{$file}}{@rownames};
}
