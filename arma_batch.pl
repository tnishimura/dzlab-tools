#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/lib";
use autodie;
use Config::General qw(ParseConfig);
use Data::Dumper;
use feature 'say';
use Getopt::Euclid qw( :vars<opt_> );
use Launch;
use Pod::Usage;
use DZUtil qw/localize common_prefix common_suffix timestamp/;
use Parallel::ForkManager;
use File::Temp qw/tempdir/;
use Cwd qw/getcwd/;
use File::Find;
use File::Basename qw/basename dirname/;
use File::Path qw/make_path remove_tree/;
use File::Spec::Functions qw/canonpath catdir catfile updir/;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) if !$opt_conf;

my $pm = Parallel::ForkManager->new($opt_parallel);

my %config = ParseConfig($opt_conf);

#my $gene = "/wip/tools/annotations/AT/gmod/TAIR8_genes.gff";
#my $transposon = "/wip/tools/annotations/AT/gmod/TAIR8_TE_mCG_min30bp_MERGED.gff";

say STDERR "building arma caches...";
while (my ($conf_name,$conf_hash) = each %config) {
    my ($binwidth,  $distance,  $stopflag,  $stopdistance,  $end,  $extractid,  $gffannotation)
    = @{$conf_hash}{'bin-width', 'distance', 'stop-flag', 'stop-distance', 'end', 'extract-id', 'gff-annotation'};

    if (! -f $gffannotation){
        die "can't find $gffannotation?"
    }

    launch("arma_build.pl -b $binwidth -d $distance -s $stopflag -k $stopdistance -g $gffannotation -x $extractid "  .
        ($end == 3 ? "-3" : $end == 5 ? "-5" : die "\$end not 3 or 5?"));
}

say STDERR "starting...";
find( 
    { 
        wanted => sub {
            # $File::Find::name - filename relative to pwd
            # $File::Find::dir  - dirname relative to pwd, eq getcwd()
            # $_                - filename relative to $File::Find::dir
            if (/^(concat-)?single-c(-\w+)?$/){
                say getcwd();
                my $armadir = "armageddon";
                my @cg;
                my @chg;
                my @chh;

                if ($opt_all){
                    @cg = grep { /[\._]all[\._]/ } glob("$_/*CG*merged*");
                    @chg = grep { /[\._]all[\._]/ } glob("$_/*CHG*merged*");
                    @chh = grep { /[\._]all[\._]/ } glob("$_/*CHH*merged*");
                }
                else{
                    @cg = grep { ! /[\._]all[\._]/ } glob("$_/*CG* $_/*cg*");
                    @chg = grep { ! /[\._]all[\._]/ } glob("$_/*CHG* $_/*chg*");
                    @chh = grep { ! /[\._]all[\._]/ } glob("$_/*CHH* $_/*chh*");
                    die "if --all, there should be exactly one of each in $_"
                    unless (@cg==1 && @chg == 1 && @chh == 1);
                }

                make_path($armadir);
                while (my ($conf_name,$conf_hash) = each %config) {

                    #$logger->info("starting run $conf_name");

                    my ($binwidth,  $distance,  $stopflag,  $stopdistance,  $end,  $extractid,  $gffannotation,  $consolidate,  $extension)
                    = @{$conf_hash}{'bin-width', 'distance', 'stop-flag', 'stop-distance', 'end', 'extract-id', 'gff-annotation', 'consolidate', 'extension' };

                    my $cg_name = catfile($armadir, $_) . "-cg.$conf_name.gff";
                    my $chg_name = catfile($armadir, $_) . "-chg.$conf_name.gff";
                    my $chh_name = catfile($armadir, $_) . "-chh.$conf_name.gff";

                    my %context = ($cg_name => \@cg, $chg_name => \@chg, $chh_name => \@chh);

                    if ($opt_dry){
                        say $cg_name;
                        say join "\n", @cg;
                        say $chg_name;
                        say join "\n", @chg;
                        say $chh_name;
                        say join "\n", @chh;
                        say "";
                    }
                    else{
                        CONTLOOP:
                        while (my ($output,$files) = each %context) {
                            next CONTLOOP if scalar(@$files) == 0;
                            $pm->start and next;
                            launch(
                                "arma.pl -b $binwidth -d $distance -s $stopflag -k $stopdistance " . 
                                "-g $gffannotation -x $extractid -o ?? -a $output.avg " .
                                ($end == 3 ? "-3" : $end == 5 ? "-5" : die "\$end not 3 or 5?") .
                                " " . 
                                join(" ", @$files),
                                expected => $output);
                            $pm->finish; 
                        }

                    }

                }
            }
        },
        follow => 1,
    }, @opt_base_dirs);

$pm->wait_all_children;

=head1 NAME

batch_ends_analysis.pl - ...

=head1 SYNOPSIS

Usage examples:

 ends_analysis_batch.pl [-b basename] [--dry|-n] [--force|-f] --conf ends.conf [--threads 4] -d basedir1 basedir2 

=head1 OPTIONS

=over

=item  -d <dir>... | --base-dirs <dir>...

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

=item  -a | --all 

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
