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

my $help;
my $conf;
my $result = GetOptions (
    "help"     => \$help,
    "conf|c=s" => \$conf,
);
pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if ($help || !$result || !$conf);  

Log::Log4perl->easy_init( { 
    level    => $DEBUG,
    #file     => ">run.log",
    layout   => '%d{HH:mm:ss} %p> (%L) %M - %m%n',
} );
my $logger = get_logger();
my $pm = Parallel::ForkManager->new(4);

my $dry = 1;
use YAML qw/LoadFile DumpFile/;
my $config = LoadFile("run.conf");

my ($left,$right,$organism,$pairs) = @{$config}{qw/left right organism pairs/};

my $outfile_l2r = "$organism-$left-to-$right.alignment";
my $outfile_r2l = "$organism-$right-to-$left.alignment";

my %globals;
while (my ($chr,$pair) = each %$pairs) {
    my ($leftfile, $rightfile) = @$pair;
    my $prefix = "$organism-$chr-$left-vs-$right";
    my $delta = "$prefix.delta";
    my $global = "$prefix.global";
    $globals{$chr} = $global;

    $logger->info("$leftfile $rightfile $prefix $delta $global");
    $pm->start and next;
    launch("nucmer -p $prefix $leftfile $rightfile",dryrun => $dry,expected => $delta);
    launch("delta-filter -g $delta > $global",dryrun => $dry,expected => $global);

    $pm->finish;
}
$pm->wait_all_children;

my %l2r = (left => $left, right => $right);
my %r2l = (right => $left, left => $right);

while (my ($chr,$file) = each %globals) {
    open my $fh, '<', $file;
    GLOBAL:
    while (defined(my $line = <$fh>)){
        chomp $line;
        my @s = split ' ', $line;
        next GLOBAL if @s != 7;
        push @{$l2r{alignment}{uc $chr}}, [ @s[0 .. 3] ];
        push @{$r2l{alignment}{uc $chr}}, [ @s[2,3,0,1] ];
    }
    close $fh;
}

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

