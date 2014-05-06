#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/lib";
use FastaReader;
use GFF::Parser;

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "meth-sites|m=s"  => \(my $methsites_file),
    "single-c|c=s"    => \(my $single_c),
    "output-file|o=s" => \(my $output_file),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || ! $single_c || ! $methsites_file||!$output_file);  
$output_file //= $single_c =~ s/gff$/complete.gff/r;

warn "reading $methsites_file";
my %strand;
{
    my $count=0;
    open my $fh, '<:crlf', $methsites_file;
    while (defined(my $line = <$fh>)){
        warn $count if $count++ % 50000 == 0;
        chomp $line;
        my ($seqid, $pos, $strand, $context) = split /\t/, $line;
        if (lc $context ne 'cg'){
            die "sorry, not non-cg support yet, this was written for assaf and his CG needs";
        }
        $strand{lc $seqid}{$pos} = $strand if lc $context eq 'cg';
    }
    close $fh;
}

my $p = GFF::Parser->new(file => $single_c, normalize => 0);

open my $output_fh, '>', $output_file;

warn "adding strand to gff";
while (defined(my $gff = $p->next)){
    my $seqid = lc $gff->sequence;
    my $pos = $gff->start;

    if (! exists $strand{$seqid}{$pos}){
        die "imposiburu";
    }
    my $strand = delete $strand{$seqid}{$pos};

    $gff->strand($strand);

    say $output_fh $gff;
}

warn "adding leftovers to gff";
for my $seqid (keys %strand) {
    for my $pos (keys %{$strand{$seqid}}) {
        my $strand = $strand{$seqid}{$pos};
        say $output_fh join "\t", $seqid, ".", "CG", $pos, $pos, ".", $strand, ".", "c=0;t=0";
    }
}

close $output_fh;

if (0 != system("sort -f -k1,1 -k4,4n -o $output_file $output_file")){
    die "couldn't sort $output_file?";
}
=head1 gff-complete-methylation.pl 

Usage examples:

 gff-complete-methylation.pl -m methsites.txt -c single-c.gff -o output.gff

=cut

