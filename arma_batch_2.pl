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

    $pm->start and next;
    launch("arma_build.pl -b $binwidth -d $distance -s $stopflag -k $stopdistance -g $gffannotation -x $extractid "  .
        ($end == 3 ? "-3" : $end == 5 ? "-5" : die "\$end not 3 or 5?"));
    $pm->finish; 
}

$pm->wait_all_children;

say STDERR "starting...";
make_path($opt_output_directory);

for my $single_c_dir (@opt_single_c_dirs) {
    say getcwd();
    my @cg;
    my @chg;
    my @chh;

    if ($opt_all){
        @cg = grep { /[\._]all[\._]/ } glob("$single_c_dir/*CG* $single_c_dir/*cg*");
        @chg = grep { /[\._]all[\._]/ } glob("$single_c_dir/*CHG* $single_c_dir/*chg*");
        @chh = grep { /[\._]all[\._]/ } glob("$single_c_dir/*CHH* $single_c_dir/*chh*");
        die "if --all, there should be exactly one of each in $single_c_dir"
        unless (@cg==1 && @chg == 1 && @chh == 1);
    }
    else{
        @cg = grep { ! /[\._]all[\._]/ } glob("$single_c_dir/*CG*gff*");
        @chg = grep { ! /[\._]all[\._]/ } glob("$single_c_dir/*CHG*gff*");
        @chh = grep { ! /[\._]all[\._]/ } glob("$single_c_dir/*CHH*gff*");
        die "if not --all, there should be more than one of each in $single_c_dir"
        unless (@cg > 1 && @chg > 1 && @chh > 1);
    }

    while (my ($conf_name,$conf_hash) = each %config) {

        my ($binwidth,  $distance,  $stopflag,  $stopdistance,  $end,  $extractid,  $gffannotation)
        = @{$conf_hash}{'bin-width', 'distance', 'stop-flag', 'stop-distance', 'end', 'extract-id', 'gff-annotation'};

        my $cg_name = catfile($opt_output_directory, $single_c_dir) . "-cg.$conf_name.gff";
        my $chg_name = catfile($opt_output_directory, $single_c_dir) . "-chg.$conf_name.gff";
        my $chh_name = catfile($opt_output_directory, $single_c_dir) . "-chh.$conf_name.gff";

        my %context = ($cg_name => \@cg, $chg_name => \@chg, $chh_name => \@chh);

        if ($opt_dry){
            say join "\n", $cg_name, @cg, $chg_name, @chg, $chh_name, @chh, "";
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
};

$pm->wait_all_children;

=head1 NAME

batch_ends_analysis.pl - ...

=head1 SYNOPSIS

Usage examples:

 ends_analysis_batch.pl [-b basename] [--dry|-n] --conf ends.conf [--threads 4] -d basedir1 basedir2 

=head1 OPTIONS

=over

=item  -o <dir> | --output-directory <dir>

output armageddon directory

=item  -d <dir>... | --single-c-dirs <dir>...

single-c directory.

=item  -c <config_file> | --conf <config_file>

=for Euclid
    config_file.type:        readable

=item --dry | -n


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
