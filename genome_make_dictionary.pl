#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use lib "$ENV{HOME}/dzlab-tools/lib";
use Launch;
use Parallel::ForkManager;
use Log::Log4perl qw/:easy/;
use Pod::Usage;
use Getopt::Long;
use File::Basename;
use File::Spec::Functions;
use List::MoreUtils qw/all/;

sub parse_delta{
    my $fh = shift;
    local $/;

    my $line = <$fh>;

    my @accum;
    # 20682761 21247247 20754282 21318768 2 2 0
    # -427283
    # 134487
    # 0
    # 20757486 20757566 4099982 4100062 0 0 0
    # 0
    # 21195575 21195644 4202480 4202549 0 0 0
    # 0
    # 21247248 21529623 21318770 21601146 3 3 0
    # -58838
    # 215882
    # -7619
    # 0
    # 21307531 21307666 8056517 8056652 0 0 0
    # 0
    # 21307535 21307697 8054529 8054691 0 0 0
    # 0
    # 21307990 21308089 8056976 8057075 0 0 0
    # 0
    # 21312830 21312923 8087201 8087294 0 0 0
    # 0
    # 21312845 21312913 8020335 8020403 0 0 0
    # 0
    # 21529629 21934924 21601152 22006449 6 6 0
    # 1070
    # -35281
    # 19137
    # -54701
    # -44
    # -285693
    # 0

    while ($line =~ 
        m/  
        (
          ^\d+ (\h \d+) {6} $ \n  # 7 numbers
          (?: 
            ^\d+$ \n              # single number lines
          )*?                    
        )
        ^0$                       # until a 0
        /gxms){

        # The four coordinates are the start and end in the reference and the
        # start and end in the query respectively. The three digits following
        # the location coordinates are the number of errors (non-identities +
        # indels), similarity errors (non-positive match scores), and stop
        # codons (does not apply to DNA alignments, will be "0"). 

        my ($start1, $end1, $start2, $end2, $numindels, undef, undef, @indel_coords)
        = split /\s/, $1; # $1 is the outer paren

        # numindels should be the count of @indel_coords, except when there's a single indel at the border,
        # which we can detect by equality of the length of the segments
        die "not the expected number of indels?: $line" 
        if ($numindels != @indel_coords && $end2 - $start2 != $end1 - $start1);

        # positive indel # = deletion from B. 
        # negative indel # = deletion from A.
        # A = XXXABCDACBDCAC$
        # B = XXXBCCDACDCAC$
        # Delta = (4, -3, 4, 0) # zero is terminator, not in @indel_coords
        # A = XXXABC.DACBDCAC$
        # B = XXX.BCCDAC.DCAC$

        for my $coord (@indel_coords) {
            if ($coord > 0){
                push @accum, [$start1, $start1 + $coord - 1, 
                              $start2, $start2 + $coord - 1, $1];
                # there was a deletion from 2, which means we skip a coord on 1
                $start1 = $start1 + $coord + 1;
                $start2 = $start2 + $coord ;
            }
            elsif ($coord < 0){
                $coord = abs($coord);

                push @accum, [$start1, $start1 + $coord - 1, 
                              $start2, $start2 + $coord - 1, $1];
                # and vice versa
                $start1 = $start1 + $coord;
                $start2 = $start2 + $coord + 1;
            }
            else{
                die "should not be a zero here";
            }
        }
        push @accum, [$start1, $end1, $start2, $end2, $1];
    }

    for my $coords (@accum) {
        my ($start1, $end1, $start2, $end2) = @$coords;
        if ($start1-$end1 != $start2-$end2){
            die Dumper $coords;
        }
        print join(",", @$coords) . "\n";
    }
    return \@accum;
}


my $help;
my $config_file;
my $dry;
my $result = GetOptions (
    "help"     => \$help,
    "conf|c=s" => \$config_file,
    "dry|n" => \$dry,
);
pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if ($help || !$result || !$config_file);  

Log::Log4perl->easy_init( { 
    level    => $DEBUG,
    #file     => ">run.log",
    layout   => '%d{HH:mm:ss} %.1p> (%L) %m%n',
} );
my $logger = get_logger();
my $pm = Parallel::ForkManager->new(4);

use YAML qw/LoadFile DumpFile/;
my $config = LoadFile($config_file);

my ($left,$right,$organism,$pairs) = @{$config}{qw/left right organism pairs/};

my $outfile_l2r = "$organism-$left-to-$right.alignment";
my $outfile_r2l = "$organism-$right-to-$left.alignment";

my %globals;
while (my ($chr,$pair) = each %$pairs) {
    my ($leftfile, $rightfile) = @$pair;
    my $prefix = basename($config_file,'.conf') . "-$organism-$chr-$left-vs-$right";
    my $delta = "$prefix.delta";
    my $global = "$prefix.global";
    $globals{$chr} = $global;

    $logger->info("$leftfile $rightfile $prefix $delta $global");
    $pm->start and next;
    launch("nucmer -p $prefix -g 0 --noextend -f $leftfile $rightfile",dryrun => $dry,expected => $delta);
    launch("delta-filter -g $delta > $global",dryrun => $dry,expected => $global);

    $pm->finish;
}
$pm->wait_all_children;

my %l2r = (left => $left, right => $right);
my %r2l = (right => $left, left => $right);

while (my ($chr,$file) = each %globals) {
    open my $fh, '<', $file;
    GLOBAL:
    for my $coords (parse_delta($fh)) {
        for my $c (@$coords) {
            my @s = @$c;
            push @{$l2r{alignment}{uc $chr}}, [ @s[0 .. 3] ];
            push @{$r2l{alignment}{uc $chr}}, [ @s[2,3,0,1] ];
        }
    }
    close $fh;
}

say Dumper \%l2r;

DumpFile($outfile_l2r, \%l2r);
DumpFile($outfile_r2l, \%r2l);


=head1 NAME

genome_make_dictionary.pl - create a dictionary to convert chromosomes.

=head1 SYNOPSIS

Usage examples:

 genome_make_dictionary.pl -c config.file

Where the config.file looks like: 

 ---
 left: 5.0            # name of the 'from' sequence
 right: 6.1           # name of the 'to' sequence
 organism: rice       # name of organism
 
 pairs:
   chr01:             # name of sequence.  
     - 5.0/chr01.con  # left file name, relative to run dir.
     - 6.1/chr01.con  # right file name, relative to run dir.
 
   chr02:             # can be repeated as necessary
     - 5.0/chr02.con
     - 6.1/chr02.con
 
   chr03: 
     - 5.0/chr03.con
     - 6.1/chr03.con
 

=cut

