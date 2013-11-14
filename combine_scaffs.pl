#!/usr/bin/env perl
use v5.14.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/lib";
use FastaReader;

my $result = GetOptions (
    "reference|r=s" => \(my $reference),
    "prefix|p=s@"   => \(my $prefix),
    "padding|pad=i" => \(my $padding_length = 0),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || ! $reference || !@$prefix);  

my $fr = FastaReader->new(file => $reference, slurp => 1);
my %regexes = map { $_ => qr/^$_/ } @$prefix;
my %groups = map { $_ => [], } @$prefix;
my @ungrouped = ();

for my $seqid ($fr->sequence_list) {
    chomp $seqid;
    my $found_pattern = 0;
    PATTERNSEARCH:
    for my $prefix (keys %regexes) {
        my $pattern = $regexes{$prefix};
        if ($seqid =~ /$pattern/){
            $found_pattern = 1;
            push @{$groups{$prefix}}, $seqid;
            last PATTERNSEARCH;
        }
    }
    if (! $found_pattern){
        push @ungrouped, $seqid;
    }
}
#######################################################################
# output

my $basename = $reference;
$basename =~ s/\.(fa|fas|fasta)$//i;

my $output_reference = IO::File->new("$basename.concat.fasta", 'w');
my $output_gff       = IO::File->new("$basename.concat.gff", 'w');

# print ungrouped chromosomes as-is
for my $s (@ungrouped) {
    print_pretty($output_reference, $s, $fr->get($s, undef, undef));
}

# print grouped chromosomes together.
my $padding = "N" x $padding_length;

for my $group_name (keys %groups) {
    my $pos = 1;
    my @accum;
    for my $id (@{$groups{$group_name}}) {
        my $seq = $fr->get($id, undef, undef);
        push @accum, $seq;
        push @accum, $padding;

        my $len = length($seq);
        $output_gff->print(join("\t", $group_name, '.', '.', $pos, $pos + $len - 1, '.', '.', '.', "ID=$id") . "\n");

        $pos += $len + length($padding);
    }
    print_pretty($output_reference, $group_name, join "", @accum);
}

$output_reference->close();
$output_gff->close();

sub print_pretty{
    my ($fh, $id, $seq) = @_;
    $fh->print(">$id\n");

    my $length = length $seq;

    my $pos = 0;
    while ($pos < $length){
        $fh->print(substr $seq, $pos, 80); 
        $fh->print("\n");
        $pos+=80;
    }
}
=head1 combine_scaffs.pl 

Combine scaffolds in a fasta file with a certain prefix.  -pad adds N's between scaffold. Creates 

 combine_scaffs.pl -r scaffolds.fasta -p prefix -pad 20

The above command creates scaffolds.concat.fasta and an "annotation" file scaffolds.concat.gff

=cut

