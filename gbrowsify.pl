#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use Launch;
use Log::Log4perl qw/get_logger/;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use DZUtil qw/timestamp/;
use File::Path qw/make_path/;
use File::Spec::Functions qw/catfile/;
use File::Basename qw/basename/;

END {close STDOUT}
$| = 1;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || ! defined $opt_reference || ! defined $opt_output_dir;

if (-e $opt_output_dir && ! -d $opt_output_dir){
    die "$opt_output_dir exists but isn't a directory?";
}
elsif (! -e $opt_output_dir){
    make_path $opt_output_dir;
}
die "can't create/write to $opt_output_dir?" if ! -d $opt_output_dir;

my $log    = catfile($opt_output_dir, join ".", "log", timestamp(), "txt");
my $sqlite = catfile($opt_output_dir, "features.sqlite");
my $refgff = catfile($opt_output_dir, basename($opt_reference, qw/fa fas fasta FA FASTA FAS/) . "gff");

my $conf=qq/
    log4perl.logger          = DEBUG, Print
    log4perl.logger.Script   = DEBUG, Print, File
    log4perl.logger.PipeLine = DEBUG, Print, File
    log4perl.logger.Module   = DEBUG, Print, File

    log4perl.appender.Print        = Log::Log4perl::Appender::Screen
    log4perl.appender.Print.layout = PatternLayout
    log4perl.appender.Print.layout.ConversionPattern = %d{HH:mm:ss} %p> %M - %m%n

    log4perl.appender.File          = Log::Log4perl::Appender::File
    log4perl.appender.File.filename = $log
    log4perl.appender.File.layout   = PatternLayout
    log4perl.appender.File.syswrite = 1
    log4perl.appender.File.layout.ConversionPattern = %d{HH:mm:ss} %p> (%L) %M - %m%n
/;
Log::Log4perl::init( \$conf );

my $logger = get_logger();

$logger->info("output dir: $opt_output_dir");
$logger->info("reference: $opt_reference");
$logger->info(defined $opt_annotation ? "annotation: $opt_annotation" : "no annotation");

launch("fasget.pl -g $opt_reference > $refgff");

my @wiggffs;
for my $methylation_file (@opt_methylation) {
    my $wiggff = catfile($opt_output_dir, basename($methylation_file, qw/gff GFF/) . "wig.gff");
    launch("methylgff2wiggle.pl -d $opt_output_dir -i $methylation_file >> $wiggff");
    push @wiggffs, $wiggff;
}
launch("bp_seqfeature_load -c -f -a DBI::SQLite -d $sqlite $opt_reference $refgff @wiggffs");

=head1 NAME

gbrowsify.pl - :

=head1 SYNOPSIS

Usage examples:

 gbrowsify.pl [options]...

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

=item --help | -h

=back

=cut
