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

#my @opt_base_dirs = split /,/, $opt_base_dirs;

$logger->info("base directories: " .  join ", ", @opt_base_dirs);

my $tempdir = tempdir(CLEANUP => 1);

for my $dir (@opt_base_dirs) {
    $logger->info("basedir  - $dir");
    if (! -d $dir){
        $logger->logdie("$dir is not a readable directory?");
    }

    # there's a possibly (from, say, a directory that was copied) that two
    # .gff.merged files will have the same name in different single-c*
    # directionries... so keep a list of ends files produced-- and check, so that
    # we can detect this situation.  Otherwise threads will operate on same files
    # at the same time and screw things up
    my %ends_file_cache; 

    for my $single_c_dir (find_single_c($dir)) {
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

            my %files;

            find( sub {
                    # $File::Find::name - filename relative to pwd
                    # $File::Find::dir  - dirname relative to pwd 
                    # $_                - filename relative to $File::Find::dir
                    if ($opt_all && /^all\.cg/i){ push @{$files{cg}}, $File::Find::name; }
                    elsif ($opt_all && /^all\.chg/i){ push @{$files{chg}}, $File::Find::name; }
                    elsif ($opt_all && /^all\.chh/i){ push @{$files{chh}}, $File::Find::name; }
                    elsif (!$opt_all && ! /_all_/ && ! /\ball\b/){ # '_all_' are the completed ones... avoid nesting
                        if    (/CG.$extension$/i)     { push @{$files{cg}}, $File::Find::name; }
                        elsif (/CHG.$extension$/i) { push @{$files{chg}}, $File::Find::name; }
                        elsif (/CHH.$extension$/i) { push @{$files{chh}}, $File::Find::name; }
                    }
                }, $single_c_dir);

            GROUPLOOP:
            while (my ($group,$group_files) = each %files) {
                $logger->info("context: $group");
                $logger->info("files:\n" . join "\n", @$group_files);
                
                my @group_bases = map {basename($_)} @$group_files;
                my $prefix = common_prefix(@group_bases);
                my $suffix = common_suffix(@group_bases);
                $prefix =~ s/\.$|\_$//;
                $suffix =~ s/^\.|\^_//;

                my $consolidated_input  = catfile($single_c_dir, $prefix) . ".all." . $suffix;
                my $ends_base = basename($consolidated_input);
                my $ends_output = catfile($ends_dir, $ends_base) . ".$conf_name.ends";
                my $avg_output  = catfile($ends_dir, $ends_base) . ".$conf_name.ends.avg";

                if (exists $ends_file_cache{$ends_output}){
                    $logger->logdie("duplicate files in single-c? script is trying to produce $ends_output twice...");
                } else{
                    $ends_file_cache{$ends_output} = 1;
                }

                # concatenate files
                if ($opt_all){
                    $consolidated_input = $group_files->[0];
                }
                elsif (! $opt_all && ! -f $consolidated_input && ! $opt_dry){
                    $logger->info("concatenating $group files into $consolidated_input");
                    open my $outfh, '>', $consolidated_input;
                    for my $f (@$group_files) {
                        open my $infh, '<', $f;
                        while (defined(my $line = <$infh>)){
                            chomp $line;
                            say $outfh $line;
                        }
                        close $infh;
                    }
                    close $outfh;
                    $logger->info("concatenation done into $consolidated_input");
                }
                else{
                    $logger->info("concatenation into $consolidated_input was already done...");
                }


                $pm->start and next GROUPLOOP;
                launch("perl -S ends_analysis.pl -g $localgff -b $binwidth -d $distance -s $stopflag -k $stopdistance "
                    .  " -x $extractid -$end -o $ends_output $consolidated_input ", 
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
}
$pm->wait_all_children;


sub find_single_c{
    my $dir = shift;
    my %accum;
    find( sub {
            if (basename($File::Find::dir) =~ /^single-c/){
                $accum{$File::Find::dir} = 1;
            }
        }, $dir);

    return keys %accum;
}


=head1 NAME

batch_ends_analysis.pl - run ends_analysis.pl in batch on bsseq or bastard outputs.

=head1 SYNOPSIS

Usage examples:

 ends_analysis_batch.pl [--dry|-n] [--force|-f] [--parallel 4] --conf ends.conf -d basedir1 basedir2 

=head1 OPTIONS

=over

=item  -d <dir>... | --base-dirs <dir>...

Root directories of bs-seq/etc run.  needs to contain a single-c* directory.

=item  -c <config_file> | --conf <config_file>

=for Euclid
    config_file.type:        readable

=item --dry | -n

=item --force | -f

=item  --all | -a

If given, assume the single-c's are pre-concatenated into files named
'all.cg*', 'all.chg*, 'all.chh*'. 

=item  --parallel <threads> | -p <threads>

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
