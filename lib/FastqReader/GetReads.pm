package FastqReader::GetReads;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use FastqReader;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(get_reads get_reads_iter);

# @HS2:90:B09PCABXX:1:1202:13269:194742 1:N:0: 
# to 
# HS2:90:B09PCABXX:1:1202:13269:194742
sub _extract_id{
    my $id_line = shift;
    $id_line =~ s/\@?([^#\s]*).*/$1/;
    return uc $id_line;
}

# given a list of read ids, return a hash of { id => sequence }
sub get_reads{
    my $file = shift;
    my @query_ids = @_;

    my $fqr = FastqReader->new(file => $file);
    return [] if ! @query_ids;

    my %lookup = map {
        _extract_id($_) => 1
    } @query_ids;

    my @results;
    while (defined(my $q = $fqr->next())){
        my ($readid, $sequence) = @$q;
        if (exists $lookup{_extract_id($readid)}){
            push @results, $q;
        }
    }
    return \@results;
}

sub get_reads_iter{
    my $file = shift;
    my @query_ids = @_;

    my %lookup = map {
        _extract_id($_) => 1
    } @query_ids;

    my $fqr = FastqReader->new(file => $file);
    return sub { return } if ! @query_ids;

    return sub {
        while (defined(my $q = $fqr->next())){
            my ($readid, $sequence) = @$q;
            if (exists $lookup{_extract_id($readid)}){
                return $q;
            }
        }
        return;
    };
}

1;

__DATA__
SINGLE:
@HS2:90:B09PCABXX:1:1202:13269:194742 1:N:0:
ATTTTTTTATTTTTGTGATCCGGTTGCGGTTTAAGTTGTTATATTTAATGATATACAGGATATTAAGTTATATTTGATTTTAAAAATTTAATTAATTTTT
+
DHHHFHHHHHHHHHFGGGGGCGGGGFFFDAGGGGGEGGGGBGGGGHDHFHHHEDHGDGDGGDGGGGEGG>GGGGG:FFFEHGHDE@GGGG@HHHHHHHHH
@HS2:90:B09PCABXX:1:1202:13292:194744 1:N:0:
GTTTAGATGTAGTGGTTGTGAAAGTTAAAATGTAGAAGGTTGATAATTTATTTGAATTAATTTGATTTTGTTGGAAGTTTGATATTTTTGGAAAAAAGTG
+
HHHHHHHHHHDGGGGGGGGGDGEDCHBHHDGGGGEDGGDBGEGGGHGHHGHHHH@HHHHHHGHHFDGGGGDDDGDEFBFFG@GEGGGDG7BDE=DB=>?@
@HS2:90:B09PCABXX:1:1202:13382:194746 1:N:0:
TGAGAGAGAAATGGATGTTGTAAAAGTTGATTGTAATTTTATTTTTGATTTGGGTTTTGTTTTAGATATTTTTTTAGGAAATTATGGTTGATTTTTGATA
+
IIIIIIEIIIHIIIIIIIIEIIIIIGFIIIIIIIIFIIIIIIIIIIIHHIIIGIEIIIHIIIIIFFIIIIIIII@@EFEFEGD@GGGGGG<GGGGGDDGD
@HS2:90:B09PCABXX:1:1202:13645:194507 1:N:0:
TTAGTTTTAATTAAAAAGTTTGGATATTTTCGGATGGTTTTAAGGGGGTTTTAAGTATTTTGGATATTTTCGATTCGGGTATTTTGGATTTTTGGATTTT
+
GGGGDGGGGEGHGDHGBGGGGGGDGHHGHHFHHHHHHFHHDDDGGBFE@DGG>GGEBGGGDDBDGGGGGGDEGGDBBDB:GGDGA@BEFFGGG83EFEEF


PAIRED:
@HWI-EAS412_0001:1:1:1031:9213#0/1
NGGCTTTTAAGATCGGGTTGCGGTTTAAGTTCTTATACTC
+HWI-EAS412_0001:1:1:1031:9213#0/1
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@HWI-EAS412_0001:1:1:1031:11200#0/1
NAACCCTAATTAGGATTTTTGAGAATGAGAGAGAGATACT
+HWI-EAS412_0001:1:1:1031:11200#0/1
KRSQSV[ZZY[Z[\]Z^^_Z\^^_BBBBBBBBBBBBBBBB
@HWI-EAS412_0001:1:1:1031:18112#0/1
NACTCAATTATACACATGACATCAAGTCATATTCGACTCT
+HWI-EAS412_0001:1:1:1031:18112#0/1
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB


