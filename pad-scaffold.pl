#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;
use FindBin;
use Pod::Usage;
use List::Util qw/first max min shuffle sum/;
use Getopt::Long;
use lib "$FindBin::Bin/lib";
use DZUtil qw/overlap/;
use GFF::Parser;
use GFF::Tree;

sub info;

my $result = GetOptions (
    "debug|d"      => \(my $debug),
    "pad-size|s=i" => \(my $pad_size = 300),
    "attribute-id|id=i" => \(my $id_tag = "Parent"),
);
my $gff_name = shift;
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || ! $gff_name);  

info("reading into tree");
my $tree = GFF::Tree->new(file => $gff_name);
my %id_to_isoforms;

info("reading again into id_to_isoforms");
{
    my $p = GFF::Parser->new(file => $gff_name);
    while (defined(my $gff = $p->next)){
        push @{$id_to_isoforms{$gff->get_attribute($id_tag)}}, $gff;
    }
}

info("now reprinting with padding");
# sort by seqid/start coord
for my $id (sort keys %id_to_isoforms) {
    my @isoforms = sort { $a->sequence cmp $b->sequence || $a->start <=> $b->start } @{$id_to_isoforms{$id}};
    info "$id has " . scalar(@isoforms) . " parts";
    die "no isoforms?" if @isoforms == 0;

    my $seqid = $isoforms[0]->sequence;
    my $isoform_start = $isoforms[0]->start;
    my $isoform_end = $isoforms[$#isoforms]->end;
    my $chromosome_start = 1;
    my $chromosome_end = 999_999_999_999; # fix me : get chromosome size from fasta and bound upstream

    # check that 
    if ($isoform_start != 1){
        my $proposed_upstream_start = max($isoform_start - $pad_size, $chromosome_start) ;
        my $proposed_upstream_end   = max($isoform_start - 1, $chromosome_start) ;

        my $accepted = 1;
        SEARCH:
        for my $result ($tree->search_overlap($seqid, $proposed_upstream_start, $proposed_upstream_end)){
            my $amount_overlap = $result->{overlap};
            my $overlapping_gff = $result->{item};
            my $overlapping_start = $overlapping_gff->{start};
            my $overlapping_end = $overlapping_gff->{end};


            if ($amount_overlap > 0){
                # check if the overlapping gff overlaps the isoform_start
                if (overlap( [$overlapping_start, $overlapping_end], [$isoform_start - 1, $isoform_start - 1])){
                    info "no acceptance because\n$overlapping_gff\noverlaps with\n$isoforms[0]\n";
                    $accepted = 0;
                    last SEARCH;
                }
                # so it must be that:
                #        proposed 
                #         [----][-------------]  isoform
                # [----------]
                # overlapping
                elsif ($proposed_upstream_start < $overlapping_end){
                    $proposed_upstream_start = $overlapping_end + 1;
                }
            }
        }
        if ($accepted){
            say join("\t", $seqid, $isoforms[0]->source // '.', 'upstream_pad', $proposed_upstream_start, $proposed_upstream_end, qw{. . .}, "$id_tag=$id",);
        }
    }
    for my $iso (@isoforms) {
        say $iso;
    }
}

#if ($debug){ say "$_ " . length($id_to_isoforms{$_}) for keys %id_to_isoforms; }

sub info{
    if ($debug){
        say STDERR @_;
    }
}
=head1 pad_scaffold.pl 

Usage examples:

 pad_scaffold.pl [options]...

=cut

