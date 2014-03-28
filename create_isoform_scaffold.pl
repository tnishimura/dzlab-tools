#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Pod::Usage;
use Getopt::Long;
use List::Util qw/first max min shuffle sum/;
use List::MoreUtils qw/all any notall uniq/;
use FindBin;
use lib "$FindBin::Bin/lib";
use FastaReader;
use GFF::Parser;

my $result = GetOptions (
    "ref|r=s" => \(my $ref_file),
    "gff|g=s" => \(my $gff_file),
    "output|o=s" => \(my $output_file),
    "locus-tag|t=s" => \(my $locus_tag = "ID"),
    "features|f=s" => \(my $features_csv),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || !$ref_file || !$gff_file || ! $output_file);  

my %features = defined $features_csv ? (map { $_ => 1} split /,/, $features_csv) : ();
my $ref = FastaReader->new(file => $ref_file, slurp => 1);
my $gffparser = GFF::Parser->new(file => $gff_file, normalize => 0);

warn "done reading ref";

open my $output_fh, '>', $output_file;
open my $annotation_fh, '>', "$output_file.gff";

# collect annotation into:
#   $isoform_to_parts{ISOFORM}{seqid}
#   $isoform_to_parts{ISOFORM}{strand}
#   $isoform_to_parts{ISOFORM}{parts} = [
#      { start => START, end => END, id => ID }
#   ]

my %isoform_to_parts;
while (defined(my $gff = $gffparser->next)){
    my $seqid      = $gff->sequence;
    my $start      = $gff->start;
    my $end        = $gff->end;
    my $strand     = $gff->strand;
    my $feature    = $gff->feature;
    my $attributes = $gff->attribute_string;

    # if feature is not on list, skip
    if (defined $feature and ! $features{$feature} ){
        next;
    }

    if ($attributes =~ /$locus_tag=( ( [^:;]+ ) :? [^;]* );?/xms){
        my $part_id = $1;
        my $isoform = $2;

        $isoform_to_parts{$isoform}{strand} //= $strand;
        if ($strand ne $isoform_to_parts{$isoform}{strand}){
            die "mutiple strands for $isoform? line $.";
        }

        $isoform_to_parts{$isoform}{seqid} //= $seqid;
        if ($seqid ne $isoform_to_parts{$isoform}{seqid}){
            die "mutiple seqid for $isoform? line $.";
        }

        push @{$isoform_to_parts{$isoform}{parts}}, {
            start => $start,
            end   => $end,
            id    => $part_id,
            feature => $feature,
        }
    }
    else{
        die "$.?";
    }
}

warn "done reading annotation into \%isoform_to_parts";

for my $isoform (keys %isoform_to_parts) {
    # warn $isoform;
    my $seqid  = $isoform_to_parts{$isoform}{seqid};
    my $strand = $isoform_to_parts{$isoform}{strand};
    my @parts  = @{$isoform_to_parts{$isoform}{parts}};
    @parts = sort { $a->{start} <=> $b->{start} } @parts;

    if ($strand eq '-'){
        @parts = reverse @parts;
    }
    my @sequence_parts = map { $ref->get($seqid, $_->{start}, $_->{end}, reverse => $strand eq '-' ? 1 : 0) } @parts;
    my $length = sum map { $_->{end} - $_->{start} + 1 } @parts;

    say $output_fh ">$isoform $length";
    say $output_fh join "", @sequence_parts;

    my $pos = 1;
    for my $p (@parts) {
        my $start = $p->{start};
        my $end = $p->{end};
        my $length = $end - $start + 1;
        my $part_id = $p->{id};
        say $annotation_fh join "\t", $isoform, ".", ".", $pos, $pos + $length - 1, ".", "+", ".", "ID=$part_id;original_start=$start;original_end=$end";
        $pos += $length;
    }
}

close $output_fh;
close $annotation_fh;

=head1 create_isoform_scaffold.pl 

This is a script like create_scaffold.pl, except that it concatenates all pieces of a gene (as determined by the locus tag) into single fasta entry.
Useful when the annotation file contatins single isoforms of each gene.
Originally written for Jessica. Example usage:

 create_isoform_scaffold.pl [-f exon,upstream_pad,downstream_pad] [-t ID] -g isoform-annotation.gff -r reference.fas -o output.fas

-f is optional, lets you define a list of features (comma separated) to grab from annotation. -t defined the locus tag (usually ID or Parent. Default ID).

=cut

__DATA__
Sample:
CHR01	MSU_osa1r7	upstream_pad	2603	2902	.	+	.	ID=LOC_Os01g01010.1:upstream_pad;Parent=LOC_Os01g01010.1
CHR01	MSU_osa1r7	exon	2903	3268	.	+	.	ID=LOC_Os01g01010.1:exon_1;Parent=LOC_Os01g01010.1
CHR01	MSU_osa1r7	exon	3354	3616	.	+	.	ID=LOC_Os01g01010.1:exon_2;Parent=LOC_Os01g01010.1
CHR01	MSU_osa1r7	exon	4357	4455	.	+	.	ID=LOC_Os01g01010.1:exon_3;Parent=LOC_Os01g01010.1
CHR01	MSU_osa1r7	exon	5457	5560	.	+	.	ID=LOC_Os01g01010.1:exon_4;Parent=LOC_Os01g01010.1
CHR01	MSU_osa1r7	exon	7136	7944	.	+	.	ID=LOC_Os01g01010.1:exon_5;Parent=LOC_Os01g01010.1
CHR01	MSU_osa1r7	exon	8028	8150	.	+	.	ID=LOC_Os01g01010.1:exon_6;Parent=LOC_Os01g01010.1
CHR01	MSU_osa1r7	exon	8232	8320	.	+	.	ID=LOC_Os01g01010.1:exon_7;Parent=LOC_Os01g01010.1
CHR01	MSU_osa1r7	exon	8408	8608	.	+	.	ID=LOC_Os01g01010.1:exon_8;Parent=LOC_Os01g01010.1
CHR01	MSU_osa1r7	exon	9210	9617	.	+	.	ID=LOC_Os01g01010.1:exon_9;Parent=LOC_Os01g01010.1
CHR01	MSU_osa1r7	exon	10104	10187	.	+	.	ID=LOC_Os01g01010.1:exon_10;Parent=LOC_Os01g01010.1
CHR01	MSU_osa1r7	exon	10274	10430	.	+	.	ID=LOC_Os01g01010.1:exon_11;Parent=LOC_Os01g01010.1
CHR01	MSU_osa1r7	exon	10504	10817	.	+	.	ID=LOC_Os01g01010.1:exon_12;Parent=LOC_Os01g01010.1
CHR01	MSU_osa1r7	downstream_pad	10818	11017	.	+	.	ID=LOC_Os01g01010.1:downstream_pad;Parent=LOC_Os01g01010.1
CHR01	MSU_osa1r7	upstream_pad	11018	11217	.	+	.	ID=LOC_Os01g01019.1:upstream_pad;Parent=LOC_Os01g01019.1
CHR01	MSU_osa1r7	exon	11218	12060	.	+	.	ID=LOC_Os01g01019.1:exon_1;Parent=LOC_Os01g01019.1
CHR01	MSU_osa1r7	exon	12152	12435	.	+	.	ID=LOC_Os01g01019.1:exon_2;Parent=LOC_Os01g01019.1
CHR01	MSU_osa1r7	downstream_pad	12436	12541	.	+	.	ID=LOC_Os01g01019.1:downstream_pad;Parent=LOC_Os01g01019.1
CHR01	MSU_osa1r7	upstream_pad	12542	12647	.	+	.	ID=LOC_Os01g01030.1:upstream_pad;Parent=LOC_Os01g01030.1
CHR01	MSU_osa1r7	exon	12648	13813	.	+	.	ID=LOC_Os01g01030.1:exon_1;Parent=LOC_Os01g01030.1
CHR01	MSU_osa1r7	exon	13906	14271	.	+	.	ID=LOC_Os01g01030.1:exon_2;Parent=LOC_Os01g01030.1
CHR01	MSU_osa1r7	exon	14359	14437	.	+	.	ID=LOC_Os01g01030.1:exon_3;Parent=LOC_Os01g01030.1
CHR01	MSU_osa1r7	exon	14969	15171	.	+	.	ID=LOC_Os01g01030.1:exon_4;Parent=LOC_Os01g01030.1
CHR01	MSU_osa1r7	exon	15266	15915	.	+	.	ID=LOC_Os01g01030.1:exon_5;Parent=LOC_Os01g01030.1
CHR01	MSU_osa1r7	downstream_pad	15916	16103	.	+	.	ID=LOC_Os01g01030.1:downstream_pad;Parent=LOC_Os01g01030.1
CHR01	MSU_osa1r7	upstream_pad	16104	16291	.	+	.	ID=LOC_Os01g01040.1:upstream_pad;Parent=LOC_Os01g01040.1
CHR01	MSU_osa1r7	exon	16292	16976	.	+	.	ID=LOC_Os01g01040.1:exon_1;Parent=LOC_Os01g01040.1
CHR01	MSU_osa1r7	exon	17383	17474	.	+	.	ID=LOC_Os01g01040.1:exon_2;Parent=LOC_Os01g01040.1
CHR01	MSU_osa1r7	exon	17558	18258	.	+	.	ID=LOC_Os01g01040.1:exon_3;Parent=LOC_Os01g01040.1
CHR01	MSU_osa1r7	exon	18501	18571	.	+	.	ID=LOC_Os01g01040.1:exon_4;Parent=LOC_Os01g01040.1
CHR01	MSU_osa1r7	exon	18968	19057	.	+	.	ID=LOC_Os01g01040.1:exon_5;Parent=LOC_Os01g01040.1
CHR01	MSU_osa1r7	exon	19142	19321	.	+	.	ID=LOC_Os01g01040.1:exon_6;Parent=LOC_Os01g01040.1
CHR01	MSU_osa1r7	exon	19531	19629	.	+	.	ID=LOC_Os01g01040.1:exon_7;Parent=LOC_Os01g01040.1
CHR01	MSU_osa1r7	exon	19734	20323	.	+	.	ID=LOC_Os01g01040.1:exon_8;Parent=LOC_Os01g01040.1
CHR01	MSU_osa1r7	downstream_pad	20324	20623	.	+	.	ID=LOC_Os01g01040.1:downstream_pad;Parent=LOC_Os01g01040.1
CHR01	MSU_osa1r7	upstream_pad	22541	22840	.	+	.	ID=LOC_Os01g01050.1:upstream_pad;Parent=LOC_Os01g01050.1
CHR01	MSU_osa1r7	exon	22841	23281	.	+	.	ID=LOC_Os01g01050.1:exon_1;Parent=LOC_Os01g01050.1
CHR01	MSU_osa1r7	exon	23572	23847	.	+	.	ID=LOC_Os01g01050.1:exon_2;Parent=LOC_Os01g01050.1
CHR01	MSU_osa1r7	exon	23962	24033	.	+	.	ID=LOC_Os01g01050.1:exon_3;Parent=LOC_Os01g01050.1
CHR01	MSU_osa1r7	exon	24492	24577	.	+	.	ID=LOC_Os01g01050.1:exon_4;Parent=LOC_Os01g01050.1
CHR01	MSU_osa1r7	exon	25445	25519	.	+	.	ID=LOC_Os01g01050.1:exon_5;Parent=LOC_Os01g01050.1
CHR01	MSU_osa1r7	exon	25883	26971	.	+	.	ID=LOC_Os01g01050.1:exon_6;Parent=LOC_Os01g01050.1
CHR01	MSU_osa1r7	downstream_pad	26972	27053	.	+	.	ID=LOC_Os01g01050.1:downstream_pad;Parent=LOC_Os01g01050.1
CHR01	MSU_osa1r7	upstream_pad	27054	27135	.	+	.	ID=LOC_Os01g01060.1:upstream_pad;Parent=LOC_Os01g01060.1
CHR01	MSU_osa1r7	exon	27136	27292	.	+	.	ID=LOC_Os01g01060.1:exon_1;Parent=LOC_Os01g01060.1
CHR01	MSU_osa1r7	exon	27370	27641	.	+	.	ID=LOC_Os01g01060.1:exon_2;Parent=LOC_Os01g01060.1
CHR01	MSU_osa1r7	exon	28090	28293	.	+	.	ID=LOC_Os01g01060.1:exon_3;Parent=LOC_Os01g01060.1
CHR01	MSU_osa1r7	exon	28365	28651	.	+	.	ID=LOC_Os01g01060.1:exon_4;Parent=LOC_Os01g01060.1
CHR01	MSU_osa1r7	downstream_pad	28652	28951	.	+	.	ID=LOC_Os01g01060.1:downstream_pad;Parent=LOC_Os01g01060.1
