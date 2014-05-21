#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Getopt::Long;
use Pod::Usage;
use IO::All;
use File::Basename qw/basename/;

use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Statistics qw/methyl_stats/;
# use YAML qw/Dump/;

my $result = GetOptions (
    "triplet|t"   => \(my $is_triplet),
    "feature|f=s" => \(my $feature),
    "output|o=s"  => \(my $output_file),
    "nuclear|n=s"      => \(my $nuclear_pattern),
    "chloroplast|c=s"  => \(my $chloroplast_pattern),
    "mitochondria|m=s" => \(my $mitochondria_pattern),
);
usage() if ! $result || ! @ARGV;

my $regexes = {
    nuclear => $nuclear_pattern      ? qr/$nuclear_pattern/i      : qr/\d+/,
    chr     => $chloroplast_pattern  ? qr/$chloroplast_pattern/i  : qr/chrc|chrpt+/i,
    mit     => $mitochondria_pattern ? qr/$mitochondria_pattern/i : qr/chrm/i,
};

my %files;

if ($is_triplet){
    if (@ARGV != 3 ){
        usage();
    }
    %files = (
        cg => $ARGV[0],
        chg => $ARGV[1],
        chh => $ARGV[2],
    );
}
elsif ($feature){
    %files = (
        $feature => $ARGV[0],
    );
    if (@ARGV != 1){
        usage();
    }
}
else{
    %files = map {
        basename($_,'.gff') => $_
    } @ARGV;
}

my ($stats, $output) = methyl_stats(%files, $regexes);

if ($output_file){
    io($output_file)->print("$output\n");
}
else{
    say $output;
}

sub usage{ pod2usage(-verbose => 2, -noperldoc => 1); }

=head1 NAME

gff-methyl-stats.pl - produce statistics about gff files

=head1 SYNOPSIS

Usage examples:

 gff-stats.pl [-d /some/temporary/directory/] file1.gff file2.gff ...

=head1 OPTIONS

=over

=item --nuclear <pattern>

Regular expression for nuclear chromosome names. default '\d+'

=item -m <pattern> | --mitochondria <pattern>

Regular expression for mitochondria chromosome name. default 'chrm'

=item -c <pattern> | --chloroplast <pattern>

Regular expression for chloroplast chromosome name. default 'chrc|chrpt'

=item  -t | --triplet

use CG/CHG/CHH as feature names. Only when 3 files given

=item  -f <f> | --feature <f>

Feature name.  Only give when file count is 1.

=item  -o <file> | --output <file>

=back

=head1 Output Columns

=over

=item context 

=item type    

=item line_count      

=item methyl_avg      

Mean C/(C+T) score per methylation site.  (Calcualte C/(C+T) for each
methylation site in the gff file, calculate the mean.)

=item c_count 

Total C in col 9.

=item t_count 

Total T in col 9.

=item coverage        

Total C + T in col 9;

=item overall_methylation     

(Total C count in GFF file)/(Total C+T count in GFF file).  This is different
from the methyl_avg scores b/c this is in total.  methyl_mean gives the same
weight to a site where c=1;t=2 as c=10;t=20, whereas methyl_total combines the
c and t's first.  

=item ct_mean 

Mean C+T score per methylation site.  (Calculate C+T for each methylation site
in the gff file, calculate the mean.)

=item ct_5%   

5% percentile of C+T score per methylation site.  (Calculate C+T for each methylation
site in the gff file, calculate 5%-tile.)

=item ct_25%  

25% percentile.

=item ct_50%  

50% percentile.

=item ct_75%  

75% percentile.

=item ct_95%

95% percentile.


=back


=cut
