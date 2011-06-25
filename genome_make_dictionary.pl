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

        my ($start1, $end1, $start2, $end2, $numindels, undef, undef, $indel_coord, @other_indels) 
        = split /\s/, $1;

        if (! all { $_ == 1 } @other_indels){
            die "more than one indel not supported yet";
        }

        if ($numindels > 0 && $indel_coord > 0){
            push @accum, [$start1, $start1+$indel_coord-1, 
            $start2, $start2+$indel_coord-1];
            push @accum, [$start1+$indel_coord+$numindels, $end1,
            $start2+$indel_coord,            $end2];
        }
        elsif ($numindels > 0 && $indel_coord < 0){
            $indel_coord *= -1;
            push @accum, [$start1, $start1+$indel_coord-1, 
            $start2, $start2+$indel_coord-1];
            push @accum, [$start1+$indel_coord,            $end1,
            $start2+$indel_coord+$numindels, $end2];
        }
        else{
            push @accum, [$start1, $end1, $start2, $end2];
        }
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

