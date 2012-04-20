#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use File::Basename qw/basename/;
use File::Path qw/make_path remove_tree/;
use File::Spec::Functions qw/catdir catfile/;
use Getopt::Long;
use Parallel::ForkManager;
use List::Util qw/sum/;
use List::MoreUtils qw/ any notall/;

use FindBin;
use lib "$ENV{HOME}/dzlab-tools/lib";
use FastqReader;
use FastaReader;
use DZUtil qw/c2t/;
use Launch qw/cast/;

END {close STDOUT}
$| = 1;

my $result = GetOptions (
    #"tdna-file|tf=s"       => \(my $tdna_file),
    #"reads-file|rf=s"      => \(my $reads_file),
    #"flank-seq-file|ff=s"  => \(my $flank_seq_file),
    "outdir|o=s"           => \(my $outdir),
    "flank-size|fs=i"      => \(my $flank_size = 30),
    "flank-mismatch|fm=i"  => \(my $flank_mismatch = 2),
    "tdna-mismatch|tm=i"   => \(my $tdna_mismatch = 2),
    "gap|g=i"              => \(my $gap = 5),
    "prefix-min-size|pm=i" => \(my $prefix_min_size = 20),
    "parallel|p=i"         => \(my $parallel = 0),
);

usage() if (!$result);

my ($tdna_file, $flank_seq_file, $reads_file) = @ARGV;

usage() if any { ! defined $_ || ! -f $_ } ($tdna_file, $flank_seq_file, $reads_file);
usage() if (! defined $outdir || ! -d $outdir);

my $pm = Parallel::ForkManager->new($parallel);

sub usage {
    say "$0 [options] -o outdir tdna-file flank-seq-file reads-file";
    exit 1;
}

#######################################################################

my $bait_file = catfile($outdir, "chromosomal-bait");
my $split_dir = catdir($outdir, "split_reads");
make_path $split_dir;
my $split_prefix = catdir($split_dir, basename($reads_file)) . ".";

#######################################################################
# create bait from flank seq 
{
    my $fr = FastaReader->new(file => $flank_seq_file);
    my @seqs_in_flank = $fr->sequence_list();
    if (! @seqs_in_flank) {
        die "flank file empty";
    }
    if (@seqs_in_flank != 1) {
        die "why are the more than one seqs in flank file?";
    }
    my $bait = substr $fr->get($seqs_in_flank[0]), 0, $flank_size;
    open my $fh, '>', $bait_file;
    say $fh ">chrbait";
    say $fh $bait;
    close $fh;
}

#######################################################################
# split the reads in smaller chunks (and convert to fasta/c2t)

my @split_reads;
SPLIT:
{
    my $max_per_file = 943_718_400; # 900MB
    my $suffix = 'aa';
    my $bytes_written = 0;;

    my @already_existing = glob("$split_prefix??");
    if (@already_existing){
        if (sum(map {-s} @already_existing) > .49 * (-s $reads_file)){
            # looks to be already exists
            say "splitting the reads file seems to be already done, skipping";
            @split_reads = @already_existing;
            last SPLIT;
        }
        else{
            say "splitting seems to have been aborted prematurely, redoing";
            unlink $_ for @already_existing;
        }
    }

    my $fqr = FastqReader->new(file => $reads_file, linesper => 4);
    open my $writer, '>', $split_prefix . $suffix;
    push @split_reads, $split_prefix . $suffix;
    while (defined(my $fq = $fqr->next())){
        my ($readid, $sequence) = @$fq;
        $sequence = c2t $sequence;

        my $to_write = ">$readid\n$sequence\n";
        print $writer $to_write;

        $bytes_written += length $to_write;

        if ($bytes_written > $max_per_file){
            close $writer;
            $suffix++;
            open $writer, '>', $split_prefix . $suffix;
            push @split_reads, $split_prefix . $suffix;
            $bytes_written = 0;
        }
    }
    close $writer;
}

#######################################################################
# bowtie-build

#say Dumper \@split_reads;

for my $split (@split_reads) {
    $pm->start and next;
    if (notall { -f && -s } map { "$split.$_" } qw/1.ebwt 2.ebwt rev.1.ebwt rev.2.ebwt/){
        cast "bowtie-build --noref $split $split";
    }
    else{
        say "bowtie-build on $split already done";
    }
    $pm->finish; 
}
$pm->wait_all_children;

exit 1;

#######################################################################

my %reads; # split_read_file => { read_id => position }

# run each bowtie in separate process, collect reads and start positions
$pm->run_on_finish(sub{ # call before calling start()
        my $ref = $_[5];
        if (defined($ref)) {  
            my ($split_read, $id_to_startpos) = @$ref;
            $reads{$split_read} = $id_to_startpos;
        } 
    });

for my $split_read (@split_reads){
    my $bowtie_output = "$bait_file-vs-" . basename($split_read) . "mm$flank_mismatch";

    $pm->start and next;
    if (! -f $bowtie_output){
        cast "bowtie --norc -f -v $flank_mismatch $split_read $bait_file > $bowtie_output",
    }

    my $return_value = [ $split_read, {} ];

    # read the 
    open my $fh, '<', $bowtie_output;
    while (defined(my $line = <$fh>)){
        chomp $line;
        my (undef, undef, $readid, $position) = split /\t/, $line;
        $return_value->[0]{$readid} = $position;
    }
    close $fh;

    $pm->finish(0, $return_value);
}
$pm->wait_all_children;

#######################################################################

while (my ($split_reads_file,$id_to_startpos) = each %reads) {
    my $fqr = FastqReader->new(file => $split_reads_file, linesper => 2);
    my $id2seq = $fqr->get_reads(keys %$id_to_startpos);
    next if (keys %$id_to_startpos == 0);

    while (my ($id,$seq) = each %$id2seq) {
        my $position = $id_to_startpos->{$id};
        if ($position - 5  > 20){
            say ">PREFIX_$id";
            say substr $seq, 0, $position;
        }
    }
}

# bowtie -f -B 1 -v 1 ref/salk-pBIN-pROK2-tdna.fa.c2t bait.prefix.fasta
