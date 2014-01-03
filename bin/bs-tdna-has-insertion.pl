#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Pod::Usage;
use Getopt::Long;
use IO::All;

use FindBin;
use lib "$FindBin::Bin/../lib";
use FastaReader;
use Run::BowtieBuild;
use FastqReader::Convert;
use Conjure;
use Sam::Alignment;

my $result = GetOptions (
    "insert-file|i=s"    => \(my $insert_file),
    "reference-file|f=s" => \(my $reference_file),
    "bootstrap|b=i"      => \(my $bootstrap = 30),
    "reads-file|r=s"     => \(my $reads_file),
    "prefix|p=s"         => \(my $prefix),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (
       ! $result 
    || ! $insert_file
    || ! $reference_file
    || ! $reads_file
    || ! $prefix
);
my $sam_file = "$prefix.sam";
my $output_file = "$prefix.counts";
my $sample_file = "$prefix.insert_plus_refsamples.fasta";

use Hash::Util qw/lock_keys unlock_keys lock_hash unlock_hash/;
use Params::Validate qw/:all/;

my @sample_seqs = create_sample(
    insert_file    => $insert_file,
    reference_file => $reference_file,
    bootstrap      => $bootstrap,
    output         => $sample_file,
);

my ($bs_sample_file) = bowtie_build(
    file    => $sample_file,
    bs      => 'c2t',
    rc      => 1,
    version => 2,
    force   => 1, # b/c the sample file is randomly generated, no cacheing
);

#######################################################################
# convert reads

my $bs_reads_file = "$reads_file.c2t";
if (! -f $bs_reads_file || -s $bs_reads_file < .25 * -s $reads_file){
    fastq_convert(to_fasta => 1, methyl => 'c2t', in => $reads_file, out => $bs_reads_file);
}

#######################################################################
# run bowtie
my %count = map { lc $_ => 0 } @sample_seqs;
lock_keys(%count);
{

    open my $samfh, '>', $sam_file;

    my $count = 0;
    conjure(
        program => "bowtie2 -f -x $bs_sample_file -U $bs_reads_file",
        on_stdout => sub {
            if (/^@/){
                say $samfh $_;
            }
            else{
                warn $count if ++$count % 50000 == 0;
                my $sam = Sam::Alignment->new($_);
                if ($sam->mapped){
                    my $seqid = $sam->seqid =~ s/^RC_//r;
                    $count{lc $seqid}++;
                    say $samfh $sam;
                }
            }
        },
        on_stderr => sub {
            say
        },
    );
    close $samfh;
}

#######################################################################
# create table

my $out = io($output_file);
while (my ($seq,$count) = each %count) {
    $out->print("$seq\t$count\n");
}

# create a fasta file with the insert, and $bootstrap random sections from the
# reference genome the same size as the insert. The goal is to align the reads
# to the sample sections + the insert to see how typical the insert is (or
# isn't).
sub create_sample{
    my %opt = validate(@_, {
            insert_file    => 1,
            reference_file => 1,
            bootstrap      => 1,
            output         => 1,
        });
    lock_keys(%opt);
    say Dumper \%opt;
    open my $fh, '>', $opt{output};

    my $insert_fr = FastaReader->new(file => $opt{insert_file}, slurp => 1,);
    my $reference_fr = FastaReader->new(file => $opt{reference_file}, slurp => 1,);

    my $insert_seq_name = $insert_fr->first_sequence();
    my $insert_length = $insert_fr->get_length($insert_seq_name);

    my %ref_sequence_lengths = $reference_fr->sequence_lengths();
    my @ref_sequences = grep { /^chr\d+$/ } keys %ref_sequence_lengths;

    say $fh ">${insert_seq_name}";
    say $fh $insert_fr->get($insert_seq_name);

    my @sample_names = ($insert_seq_name);

    for my $i (1 .. $opt{bootstrap}) {
        my $seqid = $ref_sequences[int(rand($#ref_sequences))];
        my $len = $ref_sequence_lengths{$seqid};
        my $start = int(rand($len - $insert_length - 2000)) + 1000; # don't touch 1000bp at ends 
        my $end = $start + $insert_length - 1;
        my $sample_name = "${seqid}_${start}_${end}";
        my $seq = $reference_fr->get($seqid, $start, $end);

        redo if ($seq =~ /^N*$/); # kludge - bowtie2 apparently doesn't parse fastas appropriate on extended N? 
        
        push @sample_names, $sample_name;
        say $fh ">$sample_name";
        say $fh $seq;
    }
    close $fh;
    return @sample_names;
}


=head1 bs-tdna-has-insertion.pl 

Usage examples:

 INSERT=raw/pDS-Lox-insert.fasta
 GENOME=/home/toshiro/genomes/AT/TAIR_reference.fas
 DIR=has-insert

 for i in cmt2 Jacobsen_Wt_rep2 Jacobsen_Wt_rep3 rdr2_Jacobsen
 do
    bs-tdna-has-insertion.pl -b 1000 -i $INSERT -f $GENOME -p $DIR/$i -r raw/$i.fastq
 done

=cut

