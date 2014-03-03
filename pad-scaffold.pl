#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;
use FindBin;
use Pod::Usage;
use List::Util qw/max min sum/;
use Getopt::Long;
use lib "$FindBin::Bin/lib";
use DZUtil qw/overlap/;
use GFF::Parser;
use GFF::Tree;

my $result = GetOptions (
    "debug|d"      => \(my $debug),
    "pad-size|s=i" => \(my $pad_size = 300),
    "interpad|p=i" => \(my $interpad = 0),
    "minimum-pad-size|m=i"  => \(my $minimum_pad_size = 50),
    "attribute-id|i=i" => \(my $id_tag = "Parent"),
);
my $gff_name = shift;
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || ! $gff_name);  

sub info{ if ($debug){ say STDERR @_; } }

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
my @padding;

# sort by seqid/start coord
for my $id (sort keys %id_to_isoforms) {
    my @isoforms = sort { $a->sequence cmp $b->sequence || $a->start <=> $b->start } @{$id_to_isoforms{$id}};
    info "$id has " . scalar(@isoforms) . " parts";
    die "no isoforms?" if @isoforms == 0;

    my $seqid = $isoforms[0]->sequence;
    my $strand         = $isoforms[0]->strand;
    my $is_reverse     = $strand eq '-';
    # in the code, upstream means closer to position 1 
    # in the gff, upstream means upstream from the 5' end of the isoform
    my $upstream_pad   = $is_reverse ? 'downstream_pad' : 'upstream_pad'; 
    my $downstream_pad = $is_reverse ? 'upstream_pad' : 'downstream_pad';
    my $upstream_attribute = "ID=$id:$upstream_pad;Parent=$id";
    my $downstream_attribute = "ID=$id:$downstream_pad;Parent=$id";
    my $isoform_start = $isoforms[0]->start;
    my $isoform_end = $isoforms[$#isoforms]->end;
    my $chromosome_start = 1;
    my $chromosome_end = 999_999_999_999; # fix me : get chromosome size from fasta and bound upstream

    # check that there is space upstream
    if ($isoform_start != 1){
        #my $proposed_upstream_start = max($isoform_start - $pad_size, $chromosome_start) ;
        #my $proposed_upstream_end   = max($isoform_start - 1, $chromosome_start) ;
        my $proposed_upstream_start = $isoform_start - $pad_size;
        my $proposed_upstream_end   = $isoform_start - 1;

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
            die "why is upstream flipped?" if $proposed_upstream_end < $proposed_upstream_start;
            push @padding, [ $seqid, $isoforms[0]->source, $upstream_pad, $proposed_upstream_start, $proposed_upstream_end, ".", $strand, ".", $upstream_attribute ];
        }
    }
    for my $iso (@isoforms) {
        say $iso;
    }
    if ($isoform_end != $chromosome_end){
        my $proposed_downstream_start = min($isoform_end + 1, $chromosome_end) ;
        my $proposed_downstream_end   = min($isoform_end + $pad_size, $chromosome_end) ;

        my $accepted = 1;
        SEARCH:
        for my $result ($tree->search_overlap($seqid, $proposed_downstream_start, $proposed_downstream_end)){
            my $amount_overlap = $result->{overlap};
            my $overlapping_gff = $result->{item};
            my $overlapping_start = $overlapping_gff->{start};
            my $overlapping_end = $overlapping_gff->{end};

            if ($amount_overlap > 0){
                # check if the overlapping gff overlaps the isoform_end
                if (overlap( [$overlapping_start, $overlapping_end], [$isoform_end + 1, $isoform_end + 1])){
                    info "no acceptance because\n$overlapping_gff\noverlaps with\n$isoforms[0]\n";
                    $accepted = 0;
                    last SEARCH;
                }
                # so it must be that:
                #                     proposed 
                # isoform [----------][------] 
                #                        [---------] 
                #                        overlapping

                elsif ($overlapping_start < $proposed_downstream_end){
                    $proposed_downstream_end = $overlapping_start - 1;
                }
            }
        }
        if ($accepted){
            die "why is downstream flipped?" if $proposed_downstream_end < $proposed_downstream_start;
            push @padding, [ $seqid, $isoforms[0]->source, $downstream_pad, $proposed_downstream_start, $proposed_downstream_end, ".", $strand, ".", $downstream_attribute,];
        }
    }
}

# get rid of tiny pads
@padding = grep { abs($_->[4] - $_->[3]) > 5 } @padding;

for (my $i = 0 ; $i < $#padding ; $i++){
    my $this = $padding[$i];
    my $next = $padding[$i + 1];
    # next unless ($this->[2] eq 'downstream_pad' and $next->[2] eq 'upstream_pad');
    my $this_start = $this->[3];
    my $this_end   = $this->[4];
    my $next_start = $next->[3];
    my $next_end   = $next->[4];
    if (overlap([$this_start, $this_end], [$next_start, $next_end])){
        if (
            ($this_start < $next_start && $next_end < $this_end)
            ||
            ($next_start < $this_start && $this_end < $next_end)
        ){
            die "why are there engulfings? " . join(",",@$this) . "\n" . join(",", @$next);
        }

        my $middle = int(($this_end + $next_start) / 2);
        $this->[4] = $middle - int($interpad / 2);
        $next->[3] = $middle + int($interpad / 2) + 1;
        $i++; # skip next
    }
}

for my $pad (@padding) {
    if ($pad->[4] - $pad->[3] + 1> $minimum_pad_size){
        say join "\t", @$pad;
    }
}

=head1 pad_scaffold.pl 

Add padding to exon annotation file.

 pad_scaffold.pl [--minimum-pad-size 50 | -m 50] [--pad-size 300 | -s 300] single-isoform-annotation.gff > with-padding.gff

This was written for Jessica, originally for single-isoform rice annotation files.

=cut

