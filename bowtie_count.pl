#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use FastaReader;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) if $opt_help || ! $opt_reference;

END {close STDOUT}
$| = 1;

my $fr = FastaReader->new(file => $opt_reference, ht => sub { s/^>(\S+)\s.*$/$1/}, slurp => 1);
my %counts = map { $_ => 0 } $fr->sequence_list();
my %bases = %counts;

open my $fh, '<', $opt_input;

while (defined(my $line = <$fh>)){
    chomp $line;
    my ( $read_id, $strand, $target, $coordinate, $sequence, $qualities,
        $alternatives, $snp )
    = split /\t/, $line;
    if (exists $counts{$target}){
        ++$counts{$target};
        $bases{$target}+=length $sequence;
    }
    else{
        die "$target doesn't exist in $opt_reference";
    }
}
close $fh;

if ($opt_output ne '-'){
    open my $out, '>', $opt_output;
    select $out;
}

for my $seq (sort keys %counts) {
    my $bases = $bases{$seq};
    my $count = $counts{$seq};
    my $seqlen = $fr->get_length($seq);

    say join "\t", $seq, $count, $bases, $seqlen;
}


=head1 NAME

bowtie-count.pl - count the number of reads aligning to scaffold entries.

=head1 SYNOPSIS

 bowtie-freq.pl -f scaffold.fa reads_aligned_against_that_scaffold.bowtie

=head2 inputs

Count the number of reads that align to each scaffold entry in scaffold.fa by matching the name of the scaffold entry 
against the third column of the bowtie alignment file.  For each entry in the
reference fasta, its name is taken to be the text after the '>' until the first
space.  So, if the header lines are of the form

 >atg123.1|meow 123:456:.

then that entries name is 

 atg123.1|meow

These names need to match the values in the 3rd column of the bowtie file.

=head2 output

The output has 4 columns. The first column is the scaffold entry name
('atg123.1|meow').  The second is the number of reads which mapped that entry.
The 3rd is the (number of reads aligning to the entry) * (the length of the
reads). The 4th is the length of the scaffold entry.  Unlike parse_bowtie.pl,
no normalization is done-- that's up to you.

=head1 OPTIONS

=over

=item  -f <fasta> | --reference <fasta>

=for Euclid
    fasta.type:        readable

=item <input>

=for Euclid
    input.type:        readable

=item  -o <output_file> | --output <output_file>

=for Euclid
    output_file.default:     '-'

=item --help | -h

=back

=cut
