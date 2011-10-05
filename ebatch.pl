#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/lib";
use autodie;
use Config::General qw(ParseConfig);
use Data::Dumper;
use feature 'say';
use File::Basename;
use File::Find;
use File::Spec::Functions;
use Getopt::Euclid qw( :vars<opt_> );
use Launch;
use List::Util qw/first/;
use Log::Log4perl qw/:easy/;
use Pod::Usage;
use DZUtil qw/localize common_prefix common_suffix timestamp/;
use Parallel::ForkManager;
use File::Temp qw/tempdir/;

my $pm = Parallel::ForkManager->new($opt_parallel);

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) if !$opt_conf;

my $logname = "batch_ends-" . timestamp() . ".log";

use Log::Log4perl qw/get_logger/;
my $conf=qq/
    log4perl.logger          = DEBUG, Print, File

    log4perl.appender.Print        = Log::Log4perl::Appender::Screen
    log4perl.appender.Print.layout = PatternLayout
    log4perl.appender.Print.layout.ConversionPattern = %d{HH:mm:ss} %.1p> - %m%n

    log4perl.appender.File          = Log::Log4perl::Appender::File
    log4perl.appender.File.filename = $logname
    log4perl.appender.File.layout   = PatternLayout
    log4perl.appender.File.syswrite = 1
    log4perl.appender.File.layout.ConversionPattern = %d{HH:mm:ss} %.1p> - %m%n
/;
Log::Log4perl::init( \$conf );
my $logger = get_logger();

#bin-width      = 100
#distance       = 5000
#stop-flag      = 6 
#stop-distance  = 1500
#end            = 5
#extract-id     = ID
#gff-annotation = /path/to/annotation
#consolidate    = 1
#extension      = gff.merged

#IBM_HO_2_BS/single-c/*CG.gff.merged
#IBM_HO_2_BS/single-c/IBM_HO_2_BS-chr1.single-c-CG.gff.merged  IBM_HO_2_BS/single-c/IBM_HO_2_BS-chr5.single-c-CG.gff.merged
#IBM_HO_2_BS/single-c/IBM_HO_2_BS-chr2.single-c-CG.gff.merged  IBM_HO_2_BS/single-c/IBM_HO_2_BS-chrc.single-c-CG.gff.merged
#IBM_HO_2_BS/single-c/IBM_HO_2_BS-chr3.single-c-CG.gff.merged  IBM_HO_2_BS/single-c/IBM_HO_2_BS-chrm.single-c-CG.gff.merged
#IBM_HO_2_BS/single-c/IBM_HO_2_BS-chr4.single-c-CG.gff.merged

my %config = ParseConfig($opt_conf);

$logger->info("single-c directories: " .  join ", ", @opt_single_c_dirs);

my $tempdir = tempdir(CLEANUP => 1);


for my $single_c_dir (@opt_single_c_dirs) {
    $logger->info("single-c directory - $single_c_dir");

    # create an ends/ directory along side single-c, if not already existing
    my $ends_dir = catfile(dirname($single_c_dir), "ends");
    if (! -d $ends_dir){
        mkdir $ends_dir;
        if (! -d $ends_dir){
            $logger->logdie("can't create $ends_dir?");
        }
    }

    while (my ($conf_name,$conf_hash) = each %config) {

        $logger->info("starting run $conf_name");

        my ($binwidth,  $distance,  $stopflag,  $stopdistance,  $end,  $extractid,  $gffannotation,  $consolidate,  $extension)
        = @{$conf_hash}{'bin-width', 'distance', 'stop-flag', 'stop-distance', 'end', 'extract-id', 'gff-annotation', 'consolidate', 'extension' };

        # otherwise might overwrite b/c of race.... think of something better
        my $localgff = localize($gffannotation, $tempdir); 
        my $scores = $distance * 2 / $binwidth;

        $logger->debug("\$binwidth - $binwidth");
        $logger->debug("\$distance - $distance");
        $logger->debug("\$stopflag - $stopflag");
        $logger->debug("\$stopdistance - $stopdistance");
        $logger->debug("\$end - $end");
        $logger->debug("\$extractid - $extractid");
        $logger->debug("\$gffannotation - $gffannotation");
        $logger->debug("\$consolidate - $consolidate");
        $logger->debug("\$extension - $extension");

        my @dir_contents = glob(catfile($single_c_dir,"*"));
        my ($cg) = grep { basename($_) =~ /^all\.cg/i && $_ } @dir_contents;
        my ($chg) = grep { basename($_) =~ /^all\.chg/i } @dir_contents;
        my ($chh) = grep { basename($_) =~ /^all\.chh/i } @dir_contents;

        die "can't find cg, chg, and chh concat files"
        unless 3 == grep {defined && -f} ($cg, $chg, $chh);

        my %files = (cg => $cg, chg => $chg, chh => $chh);

        GROUPLOOP:
        while (my ($group,$group_file) = each %files) {
            $logger->info("context: $group");
            $logger->info("files: $group_file");

            my $ends_base = basename($group_file);
            my $ends_output = catfile($ends_dir, $ends_base) . ".$conf_name.ends";
            my $avg_output  = catfile($ends_dir, $ends_base) . ".$conf_name.ends.avg";

            $pm->start and next GROUPLOOP;
            launch("perl -S ends_analysis.pl -g $localgff -b $binwidth -d $distance -s $stopflag -k $stopdistance "
                .  " -x $extractid -$end -o $ends_output $group_file", 
                expected => $ends_output,
                dryrun => $opt_dry,
                force => $opt_force,
            );
            launch("perl -S average_ends_new.pl -s $scores -w $binwidth -o $avg_output $ends_output", 
                expected => $avg_output,
                dryrun => $opt_dry,
                force => $opt_force,
            );
            $pm->finish;
        }
    }
}
$pm->wait_all_children;


=head1 NAME

batch_ends_analysis.pl - ...

=head1 SYNOPSIS

Usage examples:

 ends_analysis_batch.pl [-b basename] [--dry|-n] [--force|-f] --conf ends.conf [--threads 4] -d basedir1 basedir2 

=head1 OPTIONS

=over

=item  -d <dir>... | --single-c-dirs <dir>...

Root directories of bs-seq/etc run.  needs to contain a single-c* directory.

=item  -c <config_file> | --conf <config_file>

=for Euclid
    config_file.type:        readable

=item --dry | -n

=item --force | -f

=item  --parallel <threads>

Number of simultaneous ends to perform.  Default 0 for no parallelization.

=for Euclid
    threads.default:     0

=back

=head1 CONFIG

The config file can be in the following format. Each section will be run.

 <genes5>
     bin-width      = 100 
     distance       = 5000
     stop-flag      = 6 
     extract-id     = ID
     consolidate    = 1
     extension      = gff.merged
     end            = 5
     stop-distance  = 1500
     gff-annotation = http://dzlab.pmb.berkeley.edu:8080/work/annotations/AT/gmod/TAIR8_genes.gff
 </genes5>

=cut
