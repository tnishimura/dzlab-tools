#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use Launch;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use DZUtil qw/timestamp clean_basename/;
use File::Path qw/make_path/;
use File::Spec::Functions qw/catfile/;
use File::Basename qw/basename/;
use File::Temp qw/tempfile/;

END {close STDOUT}
$| = 1;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || ! defined $opt_reference || ! defined $opt_output_dir || ! defined $opt_label;

if (-e $opt_output_dir && ! -d $opt_output_dir){
    die "$opt_output_dir exists but isn't a directory?";
}
elsif (! -e $opt_output_dir){
    make_path $opt_output_dir;
}
die "can't create/write to $opt_output_dir?" if ! -d $opt_output_dir;

#######################################################################

my $log    = catfile($opt_output_dir, join ".", "log", timestamp(), "txt");
my $sqlite = catfile($opt_output_dir, "seqfeature.sqlite");
$opt_annotation //= '';
my $annonorm = $opt_annotation ? catfile($opt_output_dir, basename($opt_annotation, qw/gff GFF/) . "norm.gff") : '';
my $refnorm = catfile($opt_output_dir, basename($opt_reference, qw/fa fas fasta FA FASTA FAS/) . "norm.fasta");
my $gff = catfile($opt_output_dir, "seqfeature.gff");

#######################################################################

use Log::Log4perl qw/:easy/;
Log::Log4perl->easy_init({ level => $DEBUG, layout => '%d{HH:mm:ss} %.1p > %m%n' });
my $logger = get_logger();

$logger->info("output dir: $opt_output_dir");
$logger->info("reference: $opt_reference");
$logger->info("annotation: $opt_annotation");

#######################################################################
# normalize to lower case seq names

launch("gb-normalize-fasta.sh $opt_reference > $refnorm");
if ($opt_annotation){
    launch("gb-normalize-gff.pl $opt_annotation > $annonorm");
}

#######################################################################
# create chromosome annotation

launch("fasget.pl -g $refnorm > $gff");

#######################################################################
# convert

for my $methylation_file (@opt_methylation) {
    launch("gb-methylgff2wiggle.pl -d $opt_output_dir $methylation_file >> $gff");
}

#######################################################################
# load

launch("bp_seqfeature_load -c -f -a DBI::SQLite -d $sqlite $annonorm $refnorm $gff");

=head1 NAME

gbrowsify.pl - :

=head1 SYNOPSIS

Usage examples:

 gbrowsify.pl -o outputdir -a annotation -r genome.fasta methylation1.gff methylation2.gff

=head1 OPTIONS

=over

=item <methylation>... 

=for Euclid
    methylation.type:        readable

=item  -o <dir> | --output-dir <dir>

=item  -a <gff> | --annotation <gff>

=for Euclid
    gff.type:        readable

=item  -r <fasta> | --reference <fasta>

=for Euclid
    fasta.type:        readable

=item  -l <label> | --label <label>

=item --help | -h

=back

=cut
