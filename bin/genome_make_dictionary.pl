#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Parallel::ForkManager;
use Pod::Usage;
use Getopt::Long;
use File::Basename;
use List::MoreUtils qw/all/;
use YAML qw/LoadFile DumpFile/;
use File::Path qw/make_path/;

use FindBin;
use lib "$FindBin::Bin/../lib";
use FastaReader;

my $result = GetOptions (
    "help"     => \my $help,
    "conf|c=s" => \my $config_file,
    "tmpdir|d" => \(my $tmpdir = 'tmp'),
    "note|n=s" => \my $note,
);
pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if ($help || !$result || !$config_file);  

make_path($tmpdir);
my $config = LoadFile($config_file);

my          ($organism, $left_ecotype, $right_ecotype, $left_reference_file, $right_reference_file, $seqid_correlation) = 
@{$config}{qw/organism   left_ecotype   right_ecotype   left_reference_file   right_reference_file   seqid_correlation/};

# split reference into individual chromosomes
# seqid => filename
my %left_pieces  = split_reference($left_reference_file, $tmpdir);
my %right_pieces = split_reference($right_reference_file, $tmpdir);

my @seqid_names = sort keys %$seqid_correlation;

# prepare file names
my $conf_basename = basename($config_file, ".conf");
my %prefix = map { $_ => "$tmpdir/$_-$left_ecotype-vs-$right_ecotype" } @seqid_names;
my %global = map { $_ => "$prefix{$_}.global" } @seqid_names;
my %delta  = map { $_ => "$prefix{$_}.delta"  } @seqid_names;

# die Dumper \%prefix, \%global, \%delta;

#######################################################################
# run nucmer/delta-filter

my $pm = Parallel::ForkManager->new(4);
for my $seqid (@seqid_names) {
    my $left_piece  = $left_pieces{$seqid_correlation->{$seqid}[0]};
    my $right_piece = $right_pieces{$seqid_correlation->{$seqid}[1]};

    warn("$seqid $left_piece $right_piece $prefix{$seqid} $delta{$seqid} $global{$seqid}");

    $pm->start and next;

    run("nucmer -p $prefix{$seqid} -g 0 --noextend -f $left_piece $right_piece");
    # run("nucmer -p $prefix{$seqid} --maxmatch -g 0 --noextend -f $left_piece $right_piece");
    # run("nucmer -p $prefix{$seqid} -g 0 --noextend $left_piece $right_piece");
    run("delta-filter -g $delta{$seqid} > $global{$seqid}");

    $pm->finish;
}
$pm->wait_all_children;

my %l2r_alignment = (left => $left_ecotype, right => $right_ecotype);
my %r2l_alignment = (right => $left_ecotype, left => $right_ecotype);

# while (my ($chr,$file) = each %global) {
#     open my $fh, '<', $file;
#     say Dumper [parse_delta($fh)];
#     close $fh;
# }

while (my ($chr,$file) = each %global) {
    open my $fh, '<', $file;
    GLOBAL:
    for my $c (parse_delta($fh)) {
        my @s = @$c;
        push @{$l2r_alignment{alignment}{uc $chr}}, [ @s[0,1,2,3] ];
        push @{$r2l_alignment{alignment}{uc $chr}}, [ @s[2,3,0,1] ];
    }
    close $fh;
}

# say Dumper \%l2r_alignment;

if ($note){
    DumpFile("$organism-$left_ecotype-to-$right_ecotype.$note.alignment", \%l2r_alignment);
    DumpFile("$organism-$right_ecotype-to-$left_ecotype.$note.alignment", \%r2l_alignment);
}
else{
    DumpFile("$organism-$left_ecotype-to-$right_ecotype.alignment", \%l2r_alignment);
    DumpFile("$organism-$right_ecotype-to-$left_ecotype.alignment", \%r2l_alignment);
}

#######################################################################
#######################################################################

sub run{
    my $cmd = shift;
    say $cmd;
    system($cmd) == 0 or die "couldn't run [$cmd]";
}

sub split_reference{
    my $reference = shift;
    my $dir = shift;
    my $fr = FastaReader->new(file => $reference, slurp => 0);
    my $basename = basename($reference, '.fa', '.fasta', '.fas');
    my %accum;
    for my $seqid ($fr->sequence_list) {
        my $split = "$dir/$basename.$seqid.fasta";
        if (! -f $split){
            open my $fh, '>', $split;
            $fr->dump_pretty($fh, $seqid, $seqid);
            close $fh;
        }
        $accum{$seqid} = $split;
    }
    return %accum;
}

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
          ^\-?\d+ (\h \-?\d+) {6} $ \n  # 7 numbers
          (?: 
            ^\-?\d+$ \n              # single number lines
          )*?                    
        )
        ^0$                       # until a 0
        /gxms){

        # The four coordinates are the start and end in the reference and the
        # start and end in the query respectively. The three digits following
        # the location coordinates are the number of errors (non-identities +
        # indels), similarity errors (non-positive match scores), and stop
        # codons (does not apply to DNA alignments, will be "0"). 

        my ($start1, $end1, $start2, $end2, $numchanges, undef, undef, @indel_coords)
        = split /\s/, $1; # $1 is the outer paren

        # numchanges = number of indels + number of mutations (i think). 
        # so should be at least count of @indel_coords, except when there's a single indel at the border,
        # which we can detect by equality of the length of the segments
        die "not the expected number of indels?: $line" 
        if ($numchanges < @indel_coords && $end2 - $start2 != $end1 - $start1);

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

    # for my $coords (@accum) {
    #     my ($start1, $end1, $start2, $end2) = @$coords;
    #     if ($start1-$end1 != $start2-$end2){
    #         die Dumper $coords;
    #     }
    #     print join(",", @$coords) . "\n";
    # }
    return @accum;
}

=head1 NAME

genome_make_dictionary.pl - create a dictionary to convert chromosomes.

=head1 SYNOPSIS

Usage examples:

 genome_make_dictionary.pl -c config.file

Where the config.file looks like: 

 ---
 left_reference_file  : whatever-5.0.fasta
 right_reference_file : whatever-6.1.fasta
 left_ecotype         : 5.0
 right_ecotype        : 6.1
 
 seqid_correlation:
   chr01:             # name of sequence to use in final file
     - Chr01          # name in left reference
     - chromosome1    # name in right reference

   chr02:             
     - Chr02
     - chromosome2
 
=cut

