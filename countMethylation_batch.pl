#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use File::Spec::Functions;
use FindBin;
use lib "$FindBin::Bin/lib";
use Launch;
use Log::Log4perl qw/:easy/;
use Parallel::ForkManager;
my $pm = Parallel::ForkManager->new($opt_threads);

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if !( $opt_log_file && $opt_gff && $opt_single_c_dir && $opt_windows_dir && $opt_reference && $opt_basename );

use Log::Log4perl qw/get_logger/;
my $conf=qq/
    log4perl.logger          = DEBUG, Print, File

    log4perl.appender.Print        = Log::Log4perl::Appender::Screen
    log4perl.appender.Print.layout = PatternLayout
    log4perl.appender.Print.layout.ConversionPattern = %d{HH:mm:ss} %p> (%L) %M - %m%n

    log4perl.appender.File          = Log::Log4perl::Appender::File
    log4perl.appender.File.filename = $opt_log_file
    log4perl.appender.File.layout   = PatternLayout
    log4perl.appender.File.syswrite = 1
    log4perl.appender.File.layout.ConversionPattern = %d{HH:mm:ss} %p> (%L) %M - %m%n
/;
Log::Log4perl::init( \$conf );

my $logger = get_logger();

#######################################################################
# Split gff into chromosomes

my %split_gff;

# hackish-- split_gff made to print to split file names to stderr
# does this b/c don't chromosome names before hand.
my $split_log = "$opt_gff.splitlog";

launch("perl -S split_gff.pl --sequence all $opt_gff 2> $split_log", expected => $split_log);

open my $fh, '<', $split_log;
while (defined(my $line = <$fh>)){
    chomp $line;
    my ($chr, $file) = split /\t/, $line;
    $split_gff{$chr} = $file;
}
close $fh;

#######################################################################
# Count methylation

while (my ($chr,$split) = each %split_gff) {
    my $singlec_all          = catfile($opt_single_c_dir,$opt_basename) . ".$chr.single-c.gff";
    my %singlec_contexts = map { $_ => catfile($opt_single_c_dir,$opt_basename) . ".$chr.single-c-$_.gff" } qw/CG CHG CHH/;

    # split:            path/to/s4CxL-vs-Col.4-chr1.gff
    # singlec_all:      singlecdir/basename.chr1.single-c.gff
    # singlec_contexts: singlecdir/basename.chr1.single-c-CG.gff

    launch("perl -S countMethylation.pl --ref $opt_reference --gff $split --output $singlec_all --sort",expected => $singlec_all);
    launch("perl -S split_gff.pl --feature all $singlec_all", expected => [values %singlec_contexts]);

    while (my ($context,$singlec) = each %singlec_contexts) {
        my $merged = "$singlec.merged";
        my $window = catfile($opt_windows_dir, $opt_basename) . ".$chr.w50-$context.gff";
        $pm->start and next;
        launch("perl -S compile_gff.pl -v -o $merged $singlec", expected => $merged);
        launch("perl -S window_gff.pl -w 50 -s 50 -o $window --no-skip $merged", expected => $window);
        $pm->finish;
    }
}
$pm->wait_all_children;


=head1 NAME

countMethylation_batch.pl - ...

=head1 SYNOPSIS

Usage examples:

 countMethylation_batch.pl [options]...

=head1 OPTIONS

=over

=item  -g <gff> | --gff <gff>

=for Euclid
    gff.type:        readable

=item  -s <dir> | --single-c-dir <dir>

=for Euclid
    dir.type:        readable

=item  -w <dir> | --windows-dir <dir>

=for Euclid
    dir.type:        readable

=item  -r <ref> | --reference <ref>

=for Euclid
    ref.type:        readable

=item  -b <basename> | --basename <basename>

=item  -l <log> | --log-file <log>

=item  -t <threads> | --threads <threads>

=for Euclid
    threads.default:     0

=back

=cut

