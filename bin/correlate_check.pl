#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use FastaReader;
use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "reference|f=s" => \(my $reference),
);
pod2usage(-verbose => 1) if (!$result);  
if (!$result || ! $reference){
    say "check that correlatedPairedEnds/correlatedSingleEnds file is ok";
    say "usage: correlated_check.pl -r reference.fasta files...";
    exit 1;
}

END {close STDOUT}
$| = 1;

my $fr = FastaReader->new(file => $reference,slurp => 1);

while (defined(my $line = <ARGV>)){
    chomp $line;
    my ($seq, $read_field, $target_field,$start,$end,$strand) = (split /\t/, $line)[0,2,8,3,4,6];
    next if $seq eq '.';
    
    my $target = $target_field =~ /target=([A-Z]+)/ ? $1 : die "target malformed";
    my $read   = $read_field =~ /:([A-Z]+)$/ ? $1 : die "read malformed";

    my $original = $fr->get($seq, $start, $end, coord => $strand eq '-' ? 'r' : 'f');
    my $padded = $fr->get($seq, $start-2, $end+2, coord => $strand eq '-' ? 'r' : 'f');

    if ($padded ne $target){
        die "target mismatch: ($.)\n$line\n$padded\n$target";
    }

    $read =~ tr/cC/tT/;
    $original =~ tr/cC/tT/;

    #say join "\n", $target, "  $original", $read, $original;
    #say "";

    if ($original ne $read && numdiff($read,$original) > 3){
        die "read mismatch: ($.)\n$line\n$original\n$read";
    }
}

sub numdiff{
    my ($x,$y)=@_;
    my @x_split = split //, $x;
    my @y_split = split //, $y;
    die "len mismatch" unless @x_split == @y_split;
    my $total = 0;
    for (0..$#x_split){
        if ($x_split[$_] ne $y_split[$_]){
            $total += 1;
        }
    }
    return $total;
}

__DATA__
chr3    U/U     chr3:1147786:1147885:-:1147818#/1:TGGTAATTTTTATTGTTTTATTTGAAATTGGTTTTATTTTAAGATTAAAT    432037  432086  1       -       0       target=CGCGGTAACCTTCATTGTTTTACTTGAAATTGGTTCCATTCCAAGACTAAATGT
chr3    U/U     chr3:1147786:1147885:-:1147818#/1:GTTTGATATTATTGAATCAATATGTTGAATATTGAATATGATGATTTGAT    432087  432136  1       -       0       target=ATGTCTGATACCACTGAACCAACATGTTGAATACTGAACATGATGATTCGATCT
chr1    U/U     chr1:432951:433050:+:433046#/1:TTATAAAATATTTTTATTTTTTGGTTTATATTTTAGTTTTTTTTTTTTTT       432951  433000  1       +       0       target=TCCCACAAAATATCTTTACCCTTTGGTTTATATTTTAGTTTTTTTTTTTTTTTT
chr1    U/U     chr1:432951:433050:+:433046#/1:TTTTTTTTGTTAATTGATATAAATTAAATTTTAAAAGAGTTAATTCTGAG       433001  433050  1       +       0       target=TTTTTTTTTTGTCAACCGATATAAATTAAACCTCAAAAGAGCCAATCCCGAGGG
chr4    U/U     chr4:851894:851993:-:#/1:GTAGTAGTATTTTTGTTATGATAGTGATGGTTATATTTTTTGTTTTTTTA     727929  727978  1       -       0       target=CTGCAGTAGCATTTTCGCCATGATAGTGACGGTTATATTCTTCGTCTTCTTAAT
chr4    U/U     chr4:851894:851993:-:#/1:ATTTTTTAAGGTAATTATTATTGGAATTAAGATTTTTGGTTAAATGTTTT     727979  728028  1       -       0       target=TAATCTTTCAAGGTAATCACCATTGGAATCAAGATCTTTGGTCAAATGTTCTAT
chr1    U/U     chr1:255730:255829:+:255782#/1:AAATGTTTAATTTTGAAATTTTATGTAAATTATGTGAATATTGTGGATTG       255730  255779  1       +       0       target=ATAAATGTTCAATTTTGAAATTTTATGTAAATTATGTGAACATTGTGGATTGAA
chr1    U/U     chr1:255730:255829:+:255782#/1:AACTTATTTTAAAGTTATTAGTTATAGTTAAGTTTTATGATTAATTATTT       255780  255829  1       +       0       target=TGAACCCACTTTAAAGTTATTAGTTACAGTCAAGTCTTATGACTAATTATTTAA
chr4    U/U     chr4:1446665:1446764:+:1446668:1446674:1446715:1446729:1446749#/1:AATCTGTTGCTGTTGAGGAAATTATGGTGATAATTATTGATTGTATTTTT    1446665 1446714 1   +0       target=AAAACCTGTCGCCGTCGAGGAAATTATGGTGACAACCACCGATCGTATCTCTCT
chr4    U/U     chr4:1446665:1446764:+:1446668:1446674:1446715:1446729:1446749#/1:CTTATATTTTTGGTCTTAAGGTTTTTTTTTTTTTCATTATTTATTTTGTA    1446715 1446764 1   +0       target=CTCTTACATTCCTGGCCTCAAGGTTTTTTCTCCTCTCATCATTCATTTTGCACT

chr5	U/NM	chr5:14181068:14181167:+#/1:GTTAATAAGAGATATTAATGTGTAAATGTTTTTTTTTTTAGTATATGGAGGTTTTGGTGTTATTATTTTTTTGGAGGTTTAGAATTATATTTGAGAATTA	14181068	14181167	1	+	0	target=CTGTTAACAAGAGATACTAATGTGTAAATGTTTTTTTTCCTAGCACATGGAGGTTTTGGTGTCATTACTTTCTCGGAGGCTCAGAACTACACTTGAGAATCATC
chr5	U/NM	chr5:12508385:12508484:-#/1:AGATATGTTAATATTTATATATTAGATTTAATTATTTTTAATATGAAGATATTTTTTAAATTGGATTTTTAGTTTTATATTTTAGTAGAGAATATTTATA	14484245	14484344	1	-	0	target=ACAGATATGTCAATATCTACATATCAGATTCAATCATCTTTAACACGAAGATATTTTTCAAACCGGATTTCCAGTTTCATATCCTAGCAGAGAATATTCACAGT
chr5	U/NM	chr5:25797278:25797377:-#/2:ATGTAGTGGAGTTTTATAAGTTTTTTTTTAAGTTGATGTTAGATTTTTATTTTTTTGTTGAGTTTGTTAAAATAATAAGAGGTTTTTTTAGTGTTTAGAT	1195352	1195451	1	-	0	target=TGATGTAGCGGAGCTCCATAAGCTTTCTTCCAAGCTGACGTCAGATCCTTATCTTTTTGTTGAGTTTGTCAAAACAATAAGAGGCTTCCTCAGTGTTCAGACGG
chr4	U/NM	chr4:6874120:6874219:+#/1:GTTTAAGAATTGTAGTTAAAAATAAATTAATAATTAAATTTAAGAATTTTTGTTTGTGGGAAAAAAGATTAGGTGTTGTGTTTGAATTTTAATAAGTATA	6874120	6874219	1	+	0	target=ATGTCTAAGAATTGCAGCTAAAAATAAATTAACAATCAAATTTAAGAATTTTTGTTTGTGGGAAAAAAGACTAGGTGTCGTGTTCGAATCCTAACAAGTATAAA
chr5	U/NM	chr5:21129564:21129663:+#/2:GTTGGAGATGTAATTTTTGTAGGGAAGGATGATATAATTTTTGTTTAGAAATGAAGTTTTTTGTAATTTTATTGGTTTATGGTTTTTTAGTTAATTAAGT	21129564	21129663	1	+	0	target=TTGTTGGAGATGCAATCTCTGCAGGGAAGGACGATACAACCTTTGTCCAGAAATGAAGTTCTTTGCAACTCCACCGGTTCATGGCTCTTTAGCTAACCAAGTGG
chr1	U/NM	chr1:9572692:9572791:-#/1:TTTGTTTTATTATTAAAAATTGAAATAGTTGGAAAAAGAAAAAATAAAAGTGAGATTAAATTTAAATGGTTTGGAGAAGAAGTTAGTTTATTATTATTGA	20859773	20859872	1	-	0	target=TATTTGTTCCATCATTAAAAATCGAAATAGTCGGAAAAAGAAAAAACAAAAGCGAGACCAAACCTAAACGGTCCGGAGAAGAAGCTAGCCTATCATTACCGAAA
chr2	U/NM	chr2:9357987:9358086:+#/1:TTTGTTGGTTATGGTTTGATTGTGGTTGTTGGGTAAGATAATAAATTAGGTAGGGTTTTGTGATATTGATGTTTGTAGTTTTTAATTATGGTTAAATTGT	9357987	9358086	1	+	0	target=TATCTGTCGGTCATGGTTTGATTGTGGCTGTTGGGTAAGACAACAAATTAGGTAGGGTTTTGTGATATTGATGTCCGTAGTCCTTAATTATGGCCAAATTGTTG
chr3	U/NM	chr3:14351686:14351785:+#/1:TTTATTTTTTTTAGTAGGATTATTATTTTTTAATAAATATTTATGTGATTTGGGTGTTTTAATTAGTTTTATGTTTTTATTTGTAATTAGAAGATTGGGT	14351686	14351785	1	+	0	target=ACCCTACCTTTCTCAGTAGGACCATTATCTTTCAATAAATACCTATGTGATTTGGGTGCCTCAATTAGTCTTATGCCCCTATCTGTAACTAGAAGACTGGGTTT
chr3	U/NM	chr3:931926:932025:+#/1:TTTTAAGAGGATTTGTTGTTGATATTGAGTATGATTTTTTAAATTTTTTTAAGTTTAATATATTTGAATATAAAAGAAAATTTGTAAATATTTTTGTTTT	931926	932025	1	+	0	target=CATCTTAAGAGGATTCGTTGTTGATATTGAGTACGATTTTTCAAACTTTTCTAAGTCCAACACACTTGAATATAAAAGAAAACTTGCAAACACTCTTGTTCCAA
chr1	U/NM	chr1:12503810:12503909:-#/1:TTAGATTGGAGAATTTTATTGGTGGTGAAATTTTGTGATTGGATGATTGTATTTGTTGGATTATTGAAGGAAGATTAAATTGGGATGTTTTGGTTATTGA	17928655	17928754	1	-	0	target=GACCAGATCGGAGAATCTTACCGGCGGTGAAATTCTGCGATTGGACGATTGTATCCGTCGGATTATCGAAGGAAGACCAAACCGGGACGCTTCGGTTATTGAGG
chr1	U/NM	chr1:16011199:16011298:-#/1:GATTAGTATAGTTTGTTTTGTTGTGTTTAGTTTTTTTGTATTTGTTATATTTGTTAATTTTTATTTTTATATTTTTTGATTTTTGTTAAATAATATAAAA	14421266	14421365	1	-	0	target=TGGATTAGTACAGCTTGTCCTGTTGTGTCCAGTCTTTTTGCATCTGCCACACCTGTTAATTTTCACCTTCACACTCTTTGATTTCTGTTAAACAACATAAAATT
chr3	U/NM	chr3:23340082:23340181:-#/2:TGTTTGAAGTTTTTTGTTTTTAAAAAAAATTTAGTAATAAAATATTTGATGTTTTTTTTTAGTGATTTATTTGTTTTGTATAGAATAAAGGTAATTTTTG	130625	130724	1	-	0	target=AATGCTTGAAGCTTTTTGCCCCCAAAAAAAACTTAGTAATAAAACACTTGATGTTTTTCTTTAGCGATTCATCTGCTTTGTATAGAATAAAGGTAACCCTCGTT
chr4	U/NM	chr4:9784109:9784208:-#/1:TTTAGAAATTGGAATTTGTAGATATAAAAAAAGGAAGGATATTTATTTGAGTTTTGTTTGGGTTTTTTTTTGGTGGTTAAGTTTATAAAGTGGGTTTTAT	8800835	8800934	1	-	0	target=GCTTTAGAAATTGGAATTTGCAGACATAAAAAAAGGAAGGATATTTATTTGAGTCTTGCTTGGGCTTCTCTTTGGTGGTTAAGTTTATAAAGTGGGTCTTATAA
chr4	U/NM	chr4:3465237:3465336:-#/1:AGTTTTTTATTAGGTGGGTTGAAATTGATGTTATTTTTTATTGTGATGGTTGTGAGTTTTTAGTTGGGTAGGAGAAATTTTTTATTTATTTTTGTTTTAT	15119707	15119806	1	-	0	target=AGAGTTCTTCATTAGGCGGGTTGAAATCGATGTTACCCCTCACTGTGATGGTTGTGAGCTCTTAGTTGGGTAGGAGAAACCTCCCATCCATTTCCGCTCTATGC
chr1	U/NM	chr1:30396716:30396815:-#/1:TTTATAGAAAATTGGATGTTTAAATTGTGTAATAATTAGATATAGAGAAAGAAAAAGGTTATTATTTTATAGTGTTTGGATTTGATAAAAATATAAAATG	35749	35848	1	-	0	target=TTTTCATAGAAAATCGGATGTTCAAACTGCGTAATAATTAGATATAGAGAAAGAAAAAGGTTATTATTTTACAGTGCTTGGATTTGACAAAAATATAAAACGAC
chr1	U/NM	chr1:30090373:30090472:+#/1:TTGTTTTTTTTTTTTGTTTATATTGGATTTTTTTTTTTATTTTTTATTTGGTGTAGTTTTGTGTGTTTTAAAATTAGTGTTTTATTTATTATATTTTTTT	30090373	30090472	1	+	0	target=TGTTGCTTCTTTTTTTTGTTTATACTGGATTTTCTCTTCTATCCTCTACTTGGTGCAGCTTTGTGTGTTTTAAAACCAGTGTCTCATTCATCATACTCTTTTAC
chr3	U/NM	chr3:18359058:18359157:+#/1:TTTTTGAGAAATTTTTAGTGGGAGGAAGTGATATAAAATTAGATGATGATTTTTTTGAAGAGTAAAATTTATTTGGTGATAATATATTTTTATAAGAAGA	18359058	18359157	1	+	0	target=CCTCTCTGAGAAATTCTCAGTGGGAGGAAGTGACACAAAATCAGATGATGATTCCCTTGAAGAGTAAAATCCACTCGGTGATAACATATCTTTATAAGAAGAAT
chr3	U/NM	chr3:22405247:22405346:-#/2:GTTGTATTAGAAGTTGAGTTATGAGGAGGTTTTTGATTAGATTTAGGTTTGTATTTGAGTTGATTTGTTTTTATTTTTTTGGTTATAAATGAATTGAGTT	1065460	1065559	1	-	0	target=TCGTTGTATCAGAAGCTGAGCTATGAGGAGGCTCTTGATCAGACTCAGGCTTGTATCCGAGTTGACTCGTTTTCACCTTCTTGGTTACAAACGAATCGAGTTTC
chr3	U/NM	chr3:12028550:12028649:-#/1:GGTTTTGAGAATTTTATTGTGTTTTTTATAATGATGTGGTTATGTTTTGTAGGTAAATTTTTTTTATTTTTTATATTTTGTTATTTGTTGAGAATGGTTT	11442157	11442256	1	-	0	target=CAGGTCTTGAGAACCTCATTGTGTCTCTCACAACGATGCGGTTATGTCTTGCAGGTAAACCTTTTTCACTCCTTACACTCTGTTATTTGTTGAGAATGGTTTTT
chr5	U/NM	chr5:13847186:13847285:+#/1:AGATAAAGTATGTTGATTATTATTTAAAATTTAAATATTATTAATAAGATATAGTTAAAGTTATTAAGTGTATAGTAAAATGAAAATTTTAAGATTAAAA	13847186	13847285	1	+	0	target=ATAGACAAAGTATGTTGACTATTATTTAAAATTTAAATATCATCAATAAGATATAGTTAAAGTCATTAAGTGTATAGCAAAATGAAAATTCTAAGATTAAAATT
chr2	U/NM	chr2:14442174:14442273:-#/1:TTTAAAAAAGGGTATATTTTTTTATTTAGGAAGTAATTAGTGTATATTTGGAATGGGTTTTGAATGGTGATAATAGTTTTAAGTTGAGTAATTAGTATAG	5263087	5263186	1	-	0	target=TTTTCAAAAAAGGGTACATTCTTTTATCTAGGAAGTAACTAGTGCACACCTGGAATGGGCTTTGAATGGTGATAACAGTCCCAAGTTGAGCAACCAGCATAGAA
chr1	U/NM	chr1:16177829:16177928:+#/1:TTGAAATAAATTTTATTAAAAGTTTAGTTTTTAAATAATATTTTTTATAGATATGTTTATTATGGTATTAGTAGGATAATTTATATAAGTAATTTAAATT	16177829	16177928	1	+	0	target=GATTGAAACAAATCTTACCAAAAGTTTAGTCTTCAAACAATATCTTTCATAGATATGTTTACTACGGCATTAGTAGGATAATTTACATAAGTAATCTAAATTTA
chr5	U/NM	chr5:4242020:4242119:-#/1:AGATTAAATTTAGTTATATTTATTGAGTTGGTTTGGTGAGTTGAAATTAGGTTTATTTTTTGTGATATAATATGAAATATTTTGAATTTATTTAAAGTAG	22750610	22750709	1	-	0	target=GCAGATCAAATTCAGTTACATTCATTGAGTTGGTTTGGTGAGCCGAAATCAGGCCCATTTTCTGCGATACAACATGAAACATTCTGAACCTACCTAAAGTAGCA
chr2	U/NM	chr2:2872132:2872231:-#/2:TGAATTTGAGAGTTGTTTTTATTAATGAATTTATTTTTAGAATTTTTATTAGATTTATTTGTATTTTTTGTTGATTGTATTTTTATGTATTTGTTATTTG	16833129	16833228	1	-	0	target=GATGAACCCGAGAGTTGCCCTTATCAATGAATCCATCCTTAGAACTTCTATTAGATTTATTTGTATCCCCTGTTGATTGCACTTTTATGTACCTGTTATTTGTA
chr3	U/NM	chr3:4520296:4520395:+#/1:TTTGAGAGAGATGTTGTTTTATGTATTGTTATTATTGTTGGTTATGTTTAATTAGGTTTTGATGAAGAGGTGTTAGAGATGTTTTATAGATTGTATAGTG	4520296	4520395	1	+	0	target=TGCCTGAGAGAGATGTTGTCTCATGTACTGCTATTATTGCTGGTTATGCCCAACTAGGTCTTGATGAAGAGGCGTTAGAGATGTTTCATAGACTGCATAGTGAA
