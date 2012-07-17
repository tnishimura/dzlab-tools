#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
END {close STDOUT}
$| = 1;

my @seqid_list = qw/Chr1 Chr2 Chr3 Chr4 Chr5 Chrc Chrm/;
# my @seqid_list = qw/Chr1/;

my %tair8  = (Chr1 => 30432563, Chr2 => 19705359, Chr3 => 23470805, Chr4 => 18585042, Chr5 => 26992728, Chrc => 154478, Chrm => 366924,);
my %tair10 = (Chr1 => 30427671, Chr2 => 19698289, Chr3 => 23459830, Chr4 => 18585056, Chr5 => 26975502, Chrc => 154478, Chrm => 366924,);

my %assembly = map { $_ => [] } @seqid_list;
while (defined(my $line = <DATA>)){
    chomp $line;
    my ($seqid, $pos8, $type, $bases) = (split /\t/, $line)[0, 2, 3, 4];
    push @{$assembly{$seqid}}, [$type, $pos8, $bases];
}

my %alignment8to10 = map { $_ => [] } @seqid_list;
my %alignment10to8 = map { $_ => [] } @seqid_list;

$alignment10to8{Chrc} =$alignment8to10{Chrc} = [[(1, $tair8{Chrc}) x 2]];
$alignment10to8{Chrm} =$alignment8to10{Chrm} = [[(1, $tair8{Chrm}) x 2]];

for my $seqid (sort @seqid_list) {
    my $cursor8 = 1;
    my $cursor10 = 1;

    for (@{$assembly{$seqid}}) {
        my ($type, $pos8, $bases) = @$_;
        my $numbases = length $bases;
        if ($bases =~ /Nx(\d+)/){
            $numbases = $1;
        }
        my $length = ($pos8 - 1) - $cursor8 + 1;
        for ($type){
            when ('insertion'){
                push @{$alignment8to10{$seqid}}, [$cursor8,  $cursor8 + $length - 1,  $cursor10, $cursor10 + $length - 1];
                push @{$alignment10to8{$seqid}}, [$cursor10, $cursor10 + $length - 1, $cursor8,  $cursor8 + $length - 1];
                $cursor8  += $length;
                $cursor10 += $length + $numbases;
            }
            when ('deletion'){
                if ($cursor8 == 1){
                    $cursor8 += $numbases;
                    next;
                }
                push @{$alignment8to10{$seqid}}, [$cursor8,  $cursor8 + $length - 1,  $cursor10, $cursor10 + $length - 1];
                push @{$alignment10to8{$seqid}}, [$cursor10, $cursor10 + $length - 1, $cursor8,  $cursor8 + $length - 1];
                $cursor8  += $length + $numbases;
                $cursor10 += $length;
            }
            when ('substitution'){
                next;
            }
            default {
                die "$_???";
            }
        }
    }

    if ($cursor8 < $tair8{$seqid} && $cursor10 < $tair10{$seqid}){
        push @{$alignment8to10{$seqid}}, [$cursor8,  $tair8{$seqid},  $cursor10, $tair10{$seqid}];
        push @{$alignment10to8{$seqid}}, [$cursor10, $tair10{$seqid}, $cursor8,  $tair8{$seqid}];
    }
}

use YAML qw/Load Dump LoadFile DumpFile/;
DumpFile 'conf/at-tair8-to-tair10.alignment', {alignment =>  { map { uc $_ => $alignment8to10{$_} } keys %alignment8to10}, left => 'tair8', right => 'tair10'};
DumpFile 'conf/at-tair10-to-tair8.alignment', {alignment =>  { map { uc $_ => $alignment10to8{$_} } keys %alignment10to8}, left => 'tair10', right => 'tair8'};

=head2 Insertion @ X means:

            X
  [ ][ ][ ][ ][ ][ ][ ][ ]
  |       |\
  |       | \
  |       |  \
  [ ][ ][ ][*][ ][ ][ ][ ][ ]

=head2 Deletion @ X means:

            X
  [ ][ ][ ][*][ ][ ][ ][ ]
  |       |  /
  |       | /
  |       |/ 
  [ ][ ][ ][ ][ ][ ][ ]

=head2 Tips

* Ignore tair9 coordinates, just tair8 is enough to build dict.
* tair9->tair10 only had substituions
* subs don't affect coord

=cut

# Position on TAIR8 Assembly						***	Position on TAIR9 Assembly						
__DATA__
Chr1	TAIR8	838264	insertion	G		***	Chr1	TAIR9	838264		insertion	G
Chr1	TAIR8	859530	deletion	A		***	Chr1	TAIR9	859531	859532	deletion	A
Chr1	TAIR8	860883	insertion	A		***	Chr1	TAIR9	860883		insertion	A
Chr1	TAIR8	861073	deletion	T		***	Chr1	TAIR9	861074	861075	deletion	T
Chr1	TAIR8	869481	deletion	C		***	Chr1	TAIR9	869481	869482	deletion	C
Chr1	TAIR8	882443	insertion	T		***	Chr1	TAIR9	882442		insertion	T
Chr1	TAIR8	896233	insertion	A		***	Chr1	TAIR9	896233		insertion	A
Chr1	TAIR8	917967	deletion	G		***	Chr1	TAIR9	917968	917969	deletion	G
Chr1	TAIR8	921784	substitution	A	T	***	Chr1	TAIR9	921784		substitution	A	T
Chr1	TAIR8	936009	insertion	A		***	Chr1	TAIR9	936009		insertion	A
Chr1	TAIR8	2124250	substitution	A	G	***	Chr1	TAIR9	2124251		substitution	A	G
Chr1	TAIR8	2256361	insertion	C		***	Chr1	TAIR9	2256362		insertion	C
Chr1	TAIR8	2273166	insertion	T		***	Chr1	TAIR9	2273168		insertion	T
Chr1	TAIR8	2678388	deletion	T		***	Chr1	TAIR9	2678391	2678392	deletion	T
Chr1	TAIR8	2678582	deletion	C		***	Chr1	TAIR9	2678584	2678585	deletion	C
Chr1	TAIR8	2678594	deletion	T		***	Chr1	TAIR9	2678595	2678596	deletion	T
Chr1	TAIR8	2678700	deletion	G		***	Chr1	TAIR9	2678700	2678701	deletion	G
Chr1	TAIR8	2678735	deletion	T		***	Chr1	TAIR9	2678734	2678735	deletion	T
Chr1	TAIR8	2678751	deletion	A		***	Chr1	TAIR9	2678749	2678750	deletion	A
Chr1	TAIR8	2929300	insertion	C		***	Chr1	TAIR9	2929297		insertion	C
Chr1	TAIR8	3397482	insertion	G		***	Chr1	TAIR9	3397480		insertion	G
Chr1	TAIR8	3612993	insertion	G		***	Chr1	TAIR9	3612992		insertion	G
Chr1	TAIR8	3974329	insertion	G		***	Chr1	TAIR9	3974329		insertion	G
Chr1	TAIR8	4016007	insertion	C		***	Chr1	TAIR9	4016008		insertion	C
Chr1	TAIR8	4479486	insertion	C		***	Chr1	TAIR9	4479488		insertion	C
Chr1	TAIR8	4847916	insertion	G		***	Chr1	TAIR9	4847919		insertion	G
Chr1	TAIR8	4848549	insertion	T		***	Chr1	TAIR9	4848553		insertion	T
Chr1	TAIR8	5456138	insertion	G		***	Chr1	TAIR9	5456143		insertion	G
Chr1	TAIR8	5628144	insertion	G		***	Chr1	TAIR9	5628150		insertion	G
Chr1	TAIR8	6397361	insertion	A		***	Chr1	TAIR9	6397368		insertion	A
Chr1	TAIR8	6461647	insertion	G		***	Chr1	TAIR9	6461655		insertion	G
Chr1	TAIR8	7019232	insertion	G		***	Chr1	TAIR9	7019241		insertion	G
Chr1	TAIR8	7214815	insertion	T		***	Chr1	TAIR9	7214825		insertion	T
Chr1	TAIR8	7442693	insertion	G		***	Chr1	TAIR9	7442704		insertion	G
Chr1	TAIR8	7482221	deletion	T		***	Chr1	TAIR9	7482233	7482234	deletion	T
Chr1	TAIR8	7540565	substitution	T	C	***	Chr1	TAIR9	7540576		substitution	T	C
Chr1	TAIR8	8659326	insertion	G		***	Chr1	TAIR9	8659337		insertion	G
Chr1	TAIR8	8689012	insertion	T		***	Chr1	TAIR9	8689024		insertion	T
Chr1	TAIR8	8713017	insertion	G		***	Chr1	TAIR9	8713030		insertion	G
Chr1	TAIR8	8760282	deletion	T		***	Chr1	TAIR9	8760296	8760297	deletion	T
Chr1	TAIR8	9019335	substitution	T	A	***	Chr1	TAIR9	9019348		substitution	T	A
Chr1	TAIR8	9347806	substitution	A	G	***	Chr1	TAIR9	9347819		substitution	A	G
Chr1	TAIR8	9952606	insertion	T		***	Chr1	TAIR9	9952619		insertion	T
Chr1	TAIR8	10410815	insertion	A		***	Chr1	TAIR9	10410829		insertion	A
Chr1	TAIR8	10417364	substitution	T	A	***	Chr1	TAIR9	10417379		substitution	T	A
Chr1	TAIR8	10431782	insertion	T		***	Chr1	TAIR9	10431797		insertion	T
Chr1	TAIR8	10671380	insertion	GC		***	Chr1	TAIR9	10671396		insertion	GC
Chr1	TAIR8	10939796	insertion	T		***	Chr1	TAIR9	10939814		insertion	T
Chr1	TAIR8	11426033	insertion	A		***	Chr1	TAIR9	11426052		insertion	A
Chr1	TAIR8	12390533	insertion	C		***	Chr1	TAIR9	12390553		insertion	C
Chr1	TAIR8	12400206	insertion	C		***	Chr1	TAIR9	12400227		insertion	C
Chr1	TAIR8	12598817	substitution	A	T	***	Chr1	TAIR9	12598839		substitution	A	T
Chr1	TAIR8	12875913	substitution	C	T	***	Chr1	TAIR9	12875935		substitution	C	T
Chr1	TAIR8	13201231	deletion	Nx1234		***	Chr1	TAIR9	13201253	13201254	deletion	Nx1234
Chr1	TAIR8	13552463	substitution	A	G	***	Chr1	TAIR9	13551251		substitution	A	G
Chr1	TAIR8	13994308	deletion	Nx1229		***	Chr1	TAIR9	13993096	13993097	deletion	Nx1229
Chr1	TAIR8	14511581	substitution	C	T	***	Chr1	TAIR9	14509140		substitution	C	T
Chr1	TAIR8	15201342	substitution	A	C	***	Chr1	TAIR9	15198901		substitution	A	C
Chr1	TAIR8	15426793	substitution	G	A	***	Chr1	TAIR9	15424352		substitution	G	A
Chr1	TAIR8	15438039	substitution	C	T	***	Chr1	TAIR9	15435598		substitution	C	T
Chr1	TAIR8	15438867	substitution	A	T	***	Chr1	TAIR9	15436426		substitution	A	T
Chr1	TAIR8	15439009	substitution	A	G	***	Chr1	TAIR9	15436568		substitution	A	G
Chr1	TAIR8	15439690	substitution	C	T	***	Chr1	TAIR9	15437249		substitution	C	T
Chr1	TAIR8	15896743	insertion	C		***	Chr1	TAIR9	15894302		insertion	C
Chr1	TAIR8	16517258	substitution	C	A	***	Chr1	TAIR9	16514818		substitution	C	A
Chr1	TAIR8	16524393	substitution	T	C	***	Chr1	TAIR9	16521953		substitution	T	C
Chr1	TAIR8	17146371	substitution	C	T	***	Chr1	TAIR9	17143931		substitution	C	T
Chr1	TAIR8	17641809	deletion	Nx1229		***	Chr1	TAIR9	17639369	17639370	deletion	Nx1229
Chr1	TAIR8	17808514	insertion	G		***	Chr1	TAIR9	17804845		insertion	G
Chr1	TAIR8	18383563	substitution	G	T	***	Chr1	TAIR9	18379895		substitution	G	T
Chr1	TAIR8	19195627	deletion	A		***	Chr1	TAIR9	19191959	19191960	deletion	A
Chr1	TAIR8	19604982	insertion	G		***	Chr1	TAIR9	19601313		insertion	G
Chr1	TAIR8	20217510	insertion	C		***	Chr1	TAIR9	20213842		insertion	C
Chr1	TAIR8	20518075	substitution	T	C	***	Chr1	TAIR9	20514408		substitution	T	C
Chr1	TAIR8	20636139	insertion	G		***	Chr1	TAIR9	20632472		insertion	G
Chr1	TAIR8	20901657	insertion	G		***	Chr1	TAIR9	20897991		insertion	G
Chr1	TAIR8	21756184	substitution	A	G	***	Chr1	TAIR9	21752519		substitution	A	G
Chr1	TAIR8	23560163	insertion	CT		***	Chr1	TAIR9	23556498		insertion	CT
Chr1	TAIR8	26660141	insertion	G		***	Chr1	TAIR9	26656478		insertion	G
Chr1	TAIR8	27180129	substitution	T	C	***	Chr1	TAIR9	27176467		substitution	T	C
Chr1	TAIR8	27411422	insertion	G		***	Chr1	TAIR9	27407760		insertion	G
Chr1	TAIR8	27694265	substitution	T	C	***	Chr1	TAIR9	27690604		substitution	T	C
Chr1	TAIR8	28541151	deletion	Nx1233		***	Chr1	TAIR9	28537490	28537491	deletion	Nx1233
Chr1	TAIR8	28816720	substitution	A	G	***	Chr1	TAIR9	28811826		substitution	A	G
Chr1	TAIR8	29362270	insertion	T		***	Chr1	TAIR9	29357376		insertion	T
Chr1	TAIR8	30140448	insertion	C		***	Chr1	TAIR9	30135555		insertion	C
Chr2	TAIR8	86524	insertion	T		***	Chr2	TAIR9	86524		insertion	T
Chr2	TAIR8	275508	substitution	A	G	***	Chr2	TAIR9	275509		substitution	A	G
Chr2	TAIR8	806832	substitution	A	G	***	Chr2	TAIR9	806833		substitution	A	G
Chr2	TAIR8	937481	insertion	CC		***	Chr2	TAIR9	937482		insertion	CC
Chr2	TAIR8	973942	substitution	A	G	***	Chr2	TAIR9	973945		substitution	A	G
Chr2	TAIR8	2110641	substitution	A	G	***	Chr2	TAIR9	2110644		substitution	A	G
Chr2	TAIR8	3610369	substitution	G	T	***	Chr2	TAIR9	3610372		substitution	G	T
Chr2	TAIR8	3617842	deletion	Nx7086		***	Chr2	TAIR9	3617845	3617846	deletion	Nx7086
Chr2	TAIR8	4036906	insertion	C		***	Chr2	TAIR9	4029823		insertion	C
Chr2	TAIR8	4239858	substitution	A	G	***	Chr2	TAIR9	4232776		substitution	A	G
Chr2	TAIR8	5411744	substitution	A	G	***	Chr2	TAIR9	5404662		substitution	A	G
Chr2	TAIR8	7211571	substitution	T	C	***	Chr2	TAIR9	7204489		substitution	T	C
Chr2	TAIR8	7757998	substitution	T	C	***	Chr2	TAIR9	7750916		substitution	T	C
Chr2	TAIR8	8620442	insertion	T		***	Chr2	TAIR9	8613360		insertion	T
Chr2	TAIR8	9139759	insertion	G		***	Chr2	TAIR9	9132678		insertion	G
Chr2	TAIR8	9600266	substitution	T	C	***	Chr2	TAIR9	9593186		substitution	T	C
Chr2	TAIR8	10435327	insertion	T		***	Chr2	TAIR9	10428247		insertion	T
Chr2	TAIR8	10732510	insertion	G		***	Chr2	TAIR9	10725431		insertion	G
Chr2	TAIR8	10732513	deletion	C		***	Chr2	TAIR9	10725435	10725436	deletion	C
Chr2	TAIR8	10837163	insertion	T		***	Chr2	TAIR9	10830084		insertion	T
Chr2	TAIR8	11770298	insertion	A		***	Chr2	TAIR9	11763220		insertion	A
Chr2	TAIR8	11880270	substitution	A	T	***	Chr2	TAIR9	11873193		substitution	A	T
Chr2	TAIR8	12137157	substitution	A	G	***	Chr2	TAIR9	12130080		substitution	A	G
Chr2	TAIR8	12193566	substitution	A	G	***	Chr2	TAIR9	12186489		substitution	A	G
Chr2	TAIR8	14133729	insertion	G		***	Chr2	TAIR9	14126652		insertion	G
Chr2	TAIR8	14345167	deletion	CCC		***	Chr2	TAIR9	14338091	14338092	deletion	CCC
Chr2	TAIR8	14423736	substitution	G	T	***	Chr2	TAIR9	14416657		substitution	G	T
Chr2	TAIR8	15733660	insertion	T		***	Chr2	TAIR9	15726581		insertion	T
Chr2	TAIR8	16766678	substitution	A	G	***	Chr2	TAIR9	16759600		substitution	A	G
Chr2	TAIR8	17742280	insertion	A		***	Chr2	TAIR9	17735202		insertion	A
Chr2	TAIR8	18185997	insertion	A		***	Chr2	TAIR9	18178920		insertion	A
Chr2	TAIR8	18428218	insertion	T		***	Chr2	TAIR9	18421142		insertion	T
Chr2	TAIR8	18631524	insertion	C		***	Chr2	TAIR9	18624449		insertion	C
Chr2	TAIR8	18888957	insertion	C		***	Chr2	TAIR9	18881883		insertion	C
Chr2	TAIR8	19077693	insertion	C		***	Chr2	TAIR9	19070620		insertion	C
Chr2	TAIR8	19140468	insertion	C		***	Chr2	TAIR9	19133396		insertion	C
Chr2	TAIR8	19145750	insertion	C		***	Chr2	TAIR9	19138679		insertion	C
Chr2	TAIR8	19237565	insertion	G		***	Chr2	TAIR9	19230495		insertion	G
Chr2	TAIR8	19326107	deletion	A		***	Chr2	TAIR9	19319038	19319039	deletion	A
Chr2	TAIR8	19703781	substitution	A	G	***	Chr2	TAIR9	19696711		substitution	A	G
Chr3	TAIR8	1	deletion	Nx7		***	Chr3	TAIR9	1	2	deletion	Nx7
Chr3	TAIR8	44436	substitution	T	A	***	Chr3	TAIR9	44429		substitution	T	A
Chr3	TAIR8	1157567	insertion	C		***	Chr3	TAIR9	1157560		insertion	C
Chr3	TAIR8	2714460	insertion	A		***	Chr3	TAIR9	2714454		insertion	A
Chr3	TAIR8	2821923	deletion	T		***	Chr3	TAIR9	2821918	2821919	deletion	T
Chr3	TAIR8	2932659	insertion	T		***	Chr3	TAIR9	2932653		insertion	T
Chr3	TAIR8	3439659	deletion	T		***	Chr3	TAIR9	3439654	3439655	deletion	T
Chr3	TAIR8	3876002	deletion	T		***	Chr3	TAIR9	3875996	3875997	deletion	T
Chr3	TAIR8	4588080	substitution	T	G	***	Chr3	TAIR9	4588073		substitution	T	G
Chr3	TAIR8	5363270	insertion	G		***	Chr3	TAIR9	5363263		insertion	G
Chr3	TAIR8	5816277	substitution	C	T	***	Chr3	TAIR9	5816271		substitution	C	T
Chr3	TAIR8	7603381	deletion	C		***	Chr3	TAIR9	7603375	7603376	deletion	C
Chr3	TAIR8	9171890	deletion	Nx1230		***	Chr3	TAIR9	9171883	9171884	deletion	Nx1230
Chr3	TAIR8	9308739	substitution	A	G	***	Chr3	TAIR9	9307502		substitution	A	G
Chr3	TAIR8	9863456	substitution	A	G	***	Chr3	TAIR9	9862219		substitution	A	G
Chr3	TAIR8	10203288	substitution	G	A	***	Chr3	TAIR9	10202051		substitution	G	A
Chr3	TAIR8	10900369	substitution	A	G	***	Chr3	TAIR9	10899132		substitution	A	G
Chr3	TAIR8	11355397	deletion	Nx1242		***	Chr3	TAIR9	11354160	11354161	deletion	Nx1242
Chr3	TAIR8	13025116	deletion	Nx1234		***	Chr3	TAIR9	13022637	13022638	deletion	Nx1234
Chr3	TAIR8	13195521	substitution	G	A	***	Chr3	TAIR9	13191808		substitution	G	A
Chr3	TAIR8	13748163	deletion	Nx630		***	Chr3	TAIR9	13744450	13744451	deletion	Nx630
Chr3	TAIR8	13754065	deletion	Nx6643		***	Chr3	TAIR9	13749722	13749723	deletion	Nx6643
Chr3	TAIR8	13849978	substitution	C	T	***	Chr3	TAIR9	13838992		substitution	C	T
Chr3	TAIR8	14302875	substitution	A	G	***	Chr3	TAIR9	14291889		substitution	A	G
Chr3	TAIR8	14487543	substitution	A	T	***	Chr3	TAIR9	14476557		substitution	A	T
Chr3	TAIR8	14687919	substitution	A	G	***	Chr3	TAIR9	14676933		substitution	A	G
Chr3	TAIR8	14796015	substitution	A	G	***	Chr3	TAIR9	14785029		substitution	A	G
Chr3	TAIR8	14850783	deletion	G		***	Chr3	TAIR9	14839797	14839798	deletion	G
Chr3	TAIR8	14850823	deletion	C		***	Chr3	TAIR9	14839836	14839837	deletion	C
Chr3	TAIR8	15597340	insertion	G		***	Chr3	TAIR9	15586352		insertion	G
Chr3	TAIR8	15687459	substitution	T	C	***	Chr3	TAIR9	15676472		substitution	T	C
Chr3	TAIR8	15775867	substitution	T	C	***	Chr3	TAIR9	15764880		substitution	T	C
Chr3	TAIR8	16319839	substitution	T	C	***	Chr3	TAIR9	16308852		substitution	T	C
Chr3	TAIR8	16373573	insertion	G		***	Chr3	TAIR9	16362586		insertion	G
Chr3	TAIR8	16373574	insertion	T		***	Chr3	TAIR9	16362588		insertion	T
Chr3	TAIR8	17010818	substitution	C	A	***	Chr3	TAIR9	16999833		substitution	C	A
Chr3	TAIR8	17749194	substitution	C	T	***	Chr3	TAIR9	17738209		substitution	C	T
Chr3	TAIR8	18543471	insertion	AA		***	Chr3	TAIR9	18532486		insertion	AA
Chr3	TAIR8	18625273	insertion	G		***	Chr3	TAIR9	18614290		insertion	G
Chr3	TAIR8	18694773	insertion	C		***	Chr3	TAIR9	18683791		insertion	C
Chr3	TAIR8	18976740	insertion	T		***	Chr3	TAIR9	18965759		insertion	T
Chr3	TAIR8	19002732	insertion	A		***	Chr3	TAIR9	18991752		insertion	A
Chr3	TAIR8	19214407	insertion	C		***	Chr3	TAIR9	19203428		insertion	C
Chr3	TAIR8	19214409	insertion	C		***	Chr3	TAIR9	19203431		insertion	C
Chr3	TAIR8	19214421	insertion	C		***	Chr3	TAIR9	19203444		insertion	C
Chr3	TAIR8	19215219	deletion	T		***	Chr3	TAIR9	19204243	19204244	deletion	T
Chr3	TAIR8	19217855	substitution	A	T	***	Chr3	TAIR9	19206878		substitution	A	T
Chr3	TAIR8	19260194	substitution	T	G	***	Chr3	TAIR9	19249217		substitution	T	G
Chr3	TAIR8	19261764	deletion	C		***	Chr3	TAIR9	19250787	19250788	deletion	C
Chr3	TAIR8	20366487	insertion	G		***	Chr3	TAIR9	20355509		insertion	G
Chr3	TAIR8	20574641	deletion	T		***	Chr3	TAIR9	20563664	20563665	deletion	T
Chr3	TAIR8	20582096	deletion	C		***	Chr3	TAIR9	20571118	20571119	deletion	C
Chr3	TAIR8	20744362	substitution	A	T	***	Chr3	TAIR9	20733383		substitution	A	T
Chr3	TAIR8	21222460	insertion	A		***	Chr3	TAIR9	21211481		insertion	A
Chr3	TAIR8	21234153	substitution	T	G	***	Chr3	TAIR9	21223175		substitution	T	G
Chr3	TAIR8	21237093	insertion	G		***	Chr3	TAIR9	21226115		insertion	G
Chr3	TAIR8	21342272	substitution	A	C	***	Chr3	TAIR9	21331295		substitution	A	C
Chr3	TAIR8	22194564	insertion	C		***	Chr3	TAIR9	22183587		insertion	C
Chr3	TAIR8	22209394	insertion	G		***	Chr3	TAIR9	22198418		insertion	G
Chr3	TAIR8	23098140	deletion	G		***	Chr3	TAIR9	23087165	23087166	deletion	G
Chr3	TAIR8	23230503	insertion	C		***	Chr3	TAIR9	23219527		insertion	C
Chr4	TAIR8	1014040	substitution	C	A	***	Chr4	TAIR9	1014040		substitution	C	A
Chr4	TAIR8	1178251	substitution	T	G	***	Chr4	TAIR9	1178251		substitution	T	G
Chr4	TAIR8	1445100	insertion	C		***	Chr4	TAIR9	1445100		insertion	C
Chr4	TAIR8	1504714	insertion	C		***	Chr4	TAIR9	1504715		insertion	C
Chr4	TAIR8	2719797	insertion	A		***	Chr4	TAIR9	2719799		insertion	A
Chr4	TAIR8	3807165	substitution	T	A	***	Chr4	TAIR9	3807168		substitution	T	A
Chr4	TAIR8	5090150	substitution	A	T	***	Chr4	TAIR9	5090153		substitution	A	T
Chr4	TAIR8	6293251	insertion	C		***	Chr4	TAIR9	6293254		insertion	C
Chr4	TAIR8	6892299	insertion	A		***	Chr4	TAIR9	6892303		insertion	A
Chr4	TAIR8	6987116	insertion	C		***	Chr4	TAIR9	6987121		insertion	C
Chr4	TAIR8	7304242	substitution	A	T	***	Chr4	TAIR9	7304248		substitution	A	T
Chr4	TAIR8	7304293	substitution	G	T	***	Chr4	TAIR9	7304299		substitution	G	T
Chr4	TAIR8	7304311	substitution	A	G	***	Chr4	TAIR9	7304317		substitution	A	G
Chr4	TAIR8	7304317	substitution	A	G	***	Chr4	TAIR9	7304323		substitution	A	G
Chr4	TAIR8	7304375	substitution	C	A	***	Chr4	TAIR9	7304381		substitution	C	A
Chr4	TAIR8	7304450	deletion	G		***	Chr4	TAIR9	7304456	7304457	deletion	G
Chr4	TAIR8	7312533	deletion	A		***	Chr4	TAIR9	7312538	7312539	deletion	A
Chr4	TAIR8	7312549	deletion	C		***	Chr4	TAIR9	7312553	7312554	deletion	C
Chr4	TAIR8	7316408	deletion	G		***	Chr4	TAIR9	7316411	7316412	deletion	G
Chr4	TAIR8	7324951	insertion	G		***	Chr4	TAIR9	7324953		insertion	G
Chr4	TAIR8	7437990	substitution	T	C	***	Chr4	TAIR9	7437993		substitution	T	C
Chr4	TAIR8	7534480	deletion	G		***	Chr4	TAIR9	7534483	7534484	deletion	G
Chr4	TAIR8	7550633	insertion	A		***	Chr4	TAIR9	7550635		insertion	A
Chr4	TAIR8	7644375	insertion	C		***	Chr4	TAIR9	7644378		insertion	C
Chr4	TAIR8	8038990	insertion	T		***	Chr4	TAIR9	8038994		insertion	T
Chr4	TAIR8	8190132	substitution	A	T	***	Chr4	TAIR9	8190137		substitution	A	T
Chr4	TAIR8	8202620	substitution	C	A	***	Chr4	TAIR9	8202625		substitution	C	A
Chr4	TAIR8	8225699	deletion	C		***	Chr4	TAIR9	8225704	8225705	deletion	C
Chr4	TAIR8	8228099	substitution	T	A	***	Chr4	TAIR9	8228103		substitution	T	A
Chr4	TAIR8	8229250	substitution	T	C	***	Chr4	TAIR9	8229254		substitution	T	C
Chr4	TAIR8	8378555	insertion	A		***	Chr4	TAIR9	8378559		insertion	A
Chr4	TAIR8	8380351	deletion	T		***	Chr4	TAIR9	8380356	8380357	deletion	T
Chr4	TAIR8	8380564	insertion	T		***	Chr4	TAIR9	8380568		insertion	T
Chr4	TAIR8	8425122	substitution	A	C	***	Chr4	TAIR9	8425127		substitution	A	C
Chr4	TAIR8	8437018	substitution	T	A	***	Chr4	TAIR9	8437023		substitution	T	A
Chr4	TAIR8	8526461	deletion	C		***	Chr4	TAIR9	8526466	8526467	deletion	C
Chr4	TAIR8	8551350	deletion	G		***	Chr4	TAIR9	8551354	8551355	deletion	G
Chr4	TAIR8	8563437	insertion	G		***	Chr4	TAIR9	8563440		insertion	G
Chr4	TAIR8	8575956	insertion	A		***	Chr4	TAIR9	8575960		insertion	A
Chr4	TAIR8	8576988	insertion	A		***	Chr4	TAIR9	8576993		insertion	A
Chr4	TAIR8	8577376	insertion	G		***	Chr4	TAIR9	8577382		insertion	G
Chr4	TAIR8	8577416	deletion	G		***	Chr4	TAIR9	8577423	8577424	deletion	G
Chr4	TAIR8	8577460	deletion	A		***	Chr4	TAIR9	8577466	8577467	deletion	A
Chr4	TAIR8	8577472	substitution	T	G	***	Chr4	TAIR9	8577477		substitution	T	G
Chr4	TAIR8	8577491	substitution	C	T	***	Chr4	TAIR9	8577496		substitution	C	T
Chr4	TAIR8	8577500	deletion	T		***	Chr4	TAIR9	8577505	8577506	deletion	T
Chr4	TAIR8	8577588	substitution	A	G	***	Chr4	TAIR9	8577592		substitution	A	G
Chr4	TAIR8	8577748	insertion	A		***	Chr4	TAIR9	8577752		insertion	A
Chr4	TAIR8	8577794	insertion	A		***	Chr4	TAIR9	8577799		insertion	A
Chr4	TAIR8	8579572	deletion	T		***	Chr4	TAIR9	8579578	8579579	deletion	T
Chr4	TAIR8	8580000	substitution	A	G	***	Chr4	TAIR9	8580005		substitution	A	G
Chr4	TAIR8	8593287	deletion	C		***	Chr4	TAIR9	8593292	8593293	deletion	C
Chr4	TAIR8	8615026	deletion	T		***	Chr4	TAIR9	8615030	8615031	deletion	T
Chr4	TAIR8	8636760	substitution	T	A	***	Chr4	TAIR9	8636763		substitution	T	A
Chr4	TAIR8	8638639	deletion	C		***	Chr4	TAIR9	8638642	8638643	deletion	C
Chr4	TAIR8	8638651	insertion	T		***	Chr4	TAIR9	8638653		insertion	T
Chr4	TAIR8	8638669	deletion	T		***	Chr4	TAIR9	8638672	8638673	deletion	T
Chr4	TAIR8	8639337	deletion	A		***	Chr4	TAIR9	8639339	8639340	deletion	A
Chr4	TAIR8	8639868	deletion	G		***	Chr4	TAIR9	8639869	8639870	deletion	G
Chr4	TAIR8	8639879	deletion	G		***	Chr4	TAIR9	8639879	8639880	deletion	G
Chr4	TAIR8	8639900	substitution	G	A	***	Chr4	TAIR9	8639899		substitution	G	A
Chr4	TAIR8	8640365	substitution	G	T	***	Chr4	TAIR9	8640364		substitution	G	T
Chr4	TAIR8	8645059	substitution	G	T	***	Chr4	TAIR9	8645058		substitution	G	T
Chr4	TAIR8	8648418	substitution	G	T	***	Chr4	TAIR9	8648417		substitution	G	T
Chr4	TAIR8	8649570	deletion	T		***	Chr4	TAIR9	8649569	8649570	deletion	T
Chr4	TAIR8	8650061	insertion	G		***	Chr4	TAIR9	8650059		insertion	G
Chr4	TAIR8	8650071	insertion	G		***	Chr4	TAIR9	8650070		insertion	G
Chr4	TAIR8	8829543	insertion	G		***	Chr4	TAIR9	8829543		insertion	G
Chr4	TAIR8	8840906	insertion	A		***	Chr4	TAIR9	8840907		insertion	A
Chr4	TAIR8	8974209	deletion	T		***	Chr4	TAIR9	8974211	8974212	deletion	T
Chr4	TAIR8	8974525	deletion	A		***	Chr4	TAIR9	8974526	8974527	deletion	A
Chr4	TAIR8	8974543	substitution	G	T	***	Chr4	TAIR9	8974543		substitution	G	T
Chr4	TAIR8	8974809	substitution	A	T	***	Chr4	TAIR9	8974809		substitution	A	T
Chr4	TAIR8	8975855	deletion	G		***	Chr4	TAIR9	8975855	8975856	deletion	G
Chr4	TAIR8	8975871	deletion	T		***	Chr4	TAIR9	8975870	8975871	deletion	T
Chr4	TAIR8	8976254	deletion	T		***	Chr4	TAIR9	8976252	8976253	deletion	T
Chr4	TAIR8	8978262	insertion	C		***	Chr4	TAIR9	8978259		insertion	C
Chr4	TAIR8	8978279	insertion	C		***	Chr4	TAIR9	8978277		insertion	C
Chr4	TAIR8	8980156	deletion	C		***	Chr4	TAIR9	8980155	8980156	deletion	C
Chr4	TAIR8	8980184	deletion	A		***	Chr4	TAIR9	8980182	8980183	deletion	A
Chr4	TAIR8	8980190	deletion	A		***	Chr4	TAIR9	8980187	8980188	deletion	A
Chr4	TAIR8	8980202	deletion	C		***	Chr4	TAIR9	8980198	8980199	deletion	C
Chr4	TAIR8	8981123	substitution	C	T	***	Chr4	TAIR9	8981118		substitution	C	T
Chr4	TAIR8	8981515	substitution	G	C	***	Chr4	TAIR9	8981510		substitution	G	C
Chr4	TAIR8	8982331	deletion	C		***	Chr4	TAIR9	8982326	8982327	deletion	C
Chr4	TAIR8	8982370	deletion	A		***	Chr4	TAIR9	8982364	8982365	deletion	A
Chr4	TAIR8	8984085	deletion	A		***	Chr4	TAIR9	8984078	8984079	deletion	A
Chr4	TAIR8	8985756	substitution	G	C	***	Chr4	TAIR9	8985748		substitution	G	C
Chr4	TAIR8	8986933	deletion	C		***	Chr4	TAIR9	8986925	8986926	deletion	C
Chr4	TAIR8	8986940	deletion	A		***	Chr4	TAIR9	8986931	8986932	deletion	A
Chr4	TAIR8	8986978	deletion	CC		***	Chr4	TAIR9	8986968	8986969	deletion	CC
Chr4	TAIR8	8986997	deletion	C		***	Chr4	TAIR9	8986985	8986986	deletion	C
Chr4	TAIR8	8990753	substitution	A	C	***	Chr4	TAIR9	8990740		substitution	A	C
Chr4	TAIR8	8994912	substitution	G	T	***	Chr4	TAIR9	8994899		substitution	G	T
Chr4	TAIR8	8997483	deletion	C		***	Chr4	TAIR9	8997470	8997471	deletion	C
Chr4	TAIR8	9001460	deletion	T		***	Chr4	TAIR9	9001446	9001447	deletion	T
Chr4	TAIR8	9028713	insertion	C		***	Chr4	TAIR9	9028698		insertion	C
Chr4	TAIR8	9032486	substitution	A	T	***	Chr4	TAIR9	9032472		substitution	A	T
Chr4	TAIR8	9036427	deletion	A		***	Chr4	TAIR9	9036413	9036414	deletion	A
Chr4	TAIR8	9037055	substitution	A	G	***	Chr4	TAIR9	9037040		substitution	A	G
Chr4	TAIR8	9037148	substitution	T	G	***	Chr4	TAIR9	9037133		substitution	T	G
Chr4	TAIR8	9037215	substitution	C	A	***	Chr4	TAIR9	9037200		substitution	C	A
Chr4	TAIR8	9044641	insertion	T		***	Chr4	TAIR9	9044626		insertion	T
Chr4	TAIR8	9046804	insertion	T		***	Chr4	TAIR9	9046790		insertion	T
Chr4	TAIR8	9046824	deletion	G		***	Chr4	TAIR9	9046811	9046812	deletion	G
Chr4	TAIR8	9049395	insertion	A		***	Chr4	TAIR9	9049381		insertion	A
Chr4	TAIR8	9049396	insertion	G		***	Chr4	TAIR9	9049383		insertion	G
Chr4	TAIR8	9050451	substitution	C	A	***	Chr4	TAIR9	9050439		substitution	C	A
Chr4	TAIR8	9053620	substitution	A	C	***	Chr4	TAIR9	9053608		substitution	A	C
Chr4	TAIR8	9103911	insertion	G		***	Chr4	TAIR9	9103899		insertion	G
Chr4	TAIR8	9106127	insertion	A		***	Chr4	TAIR9	9106116		insertion	A
Chr4	TAIR8	9106604	deletion	C		***	Chr4	TAIR9	9106594	9106595	deletion	C
Chr4	TAIR8	9107806	insertion	C		***	Chr4	TAIR9	9107795		insertion	C
Chr4	TAIR8	9110772	substitution	T	A	***	Chr4	TAIR9	9110762		substitution	T	A
Chr4	TAIR8	9111354	substitution	T	A	***	Chr4	TAIR9	9111344		substitution	T	A
Chr4	TAIR8	9113228	substitution	C	T	***	Chr4	TAIR9	9113218		substitution	C	T
Chr4	TAIR8	9113674	deletion	A		***	Chr4	TAIR9	9113664	9113665	deletion	A
Chr4	TAIR8	9119323	deletion	C		***	Chr4	TAIR9	9119312	9119313	deletion	C
Chr4	TAIR8	9126055	substitution	G	T	***	Chr4	TAIR9	9126043		substitution	G	T
Chr4	TAIR8	9132741	deletion	A		***	Chr4	TAIR9	9132729	9132730	deletion	A
Chr4	TAIR8	9132759	deletion	A		***	Chr4	TAIR9	9132746	9132747	deletion	A
Chr4	TAIR8	9134812	insertion	T		***	Chr4	TAIR9	9134798		insertion	T
Chr4	TAIR8	9135715	deletion	C		***	Chr4	TAIR9	9135702	9135703	deletion	C
Chr4	TAIR8	9136925	insertion	T		***	Chr4	TAIR9	9136911		insertion	T
Chr4	TAIR8	9136931	deletion	C		***	Chr4	TAIR9	9136918	9136919	deletion	C
Chr4	TAIR8	9136959	insertion	A		***	Chr4	TAIR9	9136945		insertion	A
Chr4	TAIR8	9138766	deletion	G		***	Chr4	TAIR9	9138753	9138754	deletion	G
Chr4	TAIR8	9138814	deletion	A		***	Chr4	TAIR9	9138800	9138801	deletion	A
Chr4	TAIR8	9141479	insertion	A		***	Chr4	TAIR9	9141464		insertion	A
Chr4	TAIR8	9148836	deletion	C		***	Chr4	TAIR9	9148822	9148823	deletion	C
Chr4	TAIR8	9150641	deletion	G		***	Chr4	TAIR9	9150626	9150627	deletion	G
Chr4	TAIR8	9157229	insertion	C		***	Chr4	TAIR9	9157213		insertion	C
Chr4	TAIR8	9163898	substitution	G	T	***	Chr4	TAIR9	9163883		substitution	G	T
Chr4	TAIR8	9188055	substitution	C	T	***	Chr4	TAIR9	9188040		substitution	C	T
Chr4	TAIR8	9200408	deletion	T		***	Chr4	TAIR9	9200393	9200394	deletion	T
Chr4	TAIR8	9219309	substitution	G	C	***	Chr4	TAIR9	9219293		substitution	G	C
Chr4	TAIR8	9220035	substitution	C	T	***	Chr4	TAIR9	9220019		substitution	C	T
Chr4	TAIR8	9223300	deletion	C		***	Chr4	TAIR9	9223284	9223285	deletion	C
Chr4	TAIR8	9224064	substitution	C	T	***	Chr4	TAIR9	9224047		substitution	C	T
Chr4	TAIR8	9226778	substitution	G	C	***	Chr4	TAIR9	9226761		substitution	G	C
Chr4	TAIR8	9226804	substitution	T	A	***	Chr4	TAIR9	9226787		substitution	T	A
Chr4	TAIR8	9226835	substitution	T	A	***	Chr4	TAIR9	9226818		substitution	T	A
Chr4	TAIR8	9226838	substitution	T	A	***	Chr4	TAIR9	9226821		substitution	T	A
Chr4	TAIR8	9227318	insertion	C		***	Chr4	TAIR9	9227301		insertion	C
Chr4	TAIR8	9227388	substitution	T	G	***	Chr4	TAIR9	9227372		substitution	T	G
Chr4	TAIR8	9227418	deletion	T		***	Chr4	TAIR9	9227402	9227403	deletion	T
Chr4	TAIR8	9227427	deletion	T		***	Chr4	TAIR9	9227410	9227411	deletion	T
Chr4	TAIR8	9227443	insertion	C		***	Chr4	TAIR9	9227425		insertion	C
Chr4	TAIR8	9230475	substitution	G	A	***	Chr4	TAIR9	9230458		substitution	G	A
Chr4	TAIR8	9232541	substitution	T	A	***	Chr4	TAIR9	9232524		substitution	T	A
Chr4	TAIR8	9235202	insertion	G		***	Chr4	TAIR9	9235185		insertion	G
Chr4	TAIR8	9236210	deletion	A		***	Chr4	TAIR9	9236194	9236195	deletion	A
Chr4	TAIR8	9238839	substitution	A	G	***	Chr4	TAIR9	9238822		substitution	A	G
Chr4	TAIR8	9238840	substitution	C	A	***	Chr4	TAIR9	9238823		substitution	C	A
Chr4	TAIR8	9240376	substitution	G	A	***	Chr4	TAIR9	9240359		substitution	G	A
Chr4	TAIR8	9240491	substitution	C	A	***	Chr4	TAIR9	9240474		substitution	C	A
Chr4	TAIR8	9241274	substitution	G	A	***	Chr4	TAIR9	9241257		substitution	G	A
Chr4	TAIR8	9245644	deletion	A		***	Chr4	TAIR9	9245627	9245628	deletion	A
Chr4	TAIR8	9245669	deletion	C		***	Chr4	TAIR9	9245651	9245652	deletion	C
Chr4	TAIR8	9246479	substitution	C	T	***	Chr4	TAIR9	9246460		substitution	C	T
Chr4	TAIR8	9247899	deletion	G		***	Chr4	TAIR9	9247880	9247881	deletion	G
Chr4	TAIR8	9247932	substitution	A	G	***	Chr4	TAIR9	9247912		substitution	A	G
Chr4	TAIR8	9249528	substitution	T	C	***	Chr4	TAIR9	9249508		substitution	T	C
Chr4	TAIR8	9249866	deletion	G		***	Chr4	TAIR9	9249846	9249847	deletion	G
Chr4	TAIR8	9252090	substitution	C	A	***	Chr4	TAIR9	9252069		substitution	C	A
Chr4	TAIR8	9254678	substitution	G	A	***	Chr4	TAIR9	9254657		substitution	G	A
Chr4	TAIR8	9254734	substitution	T	A	***	Chr4	TAIR9	9254713		substitution	T	A
Chr4	TAIR8	9254738	substitution	T	G	***	Chr4	TAIR9	9254717		substitution	T	G
Chr4	TAIR8	9255164	substitution	C	G	***	Chr4	TAIR9	9255143		substitution	C	G
Chr4	TAIR8	9255205	deletion	C		***	Chr4	TAIR9	9255184	9255185	deletion	C
Chr4	TAIR8	9256739	substitution	T	A	***	Chr4	TAIR9	9256717		substitution	T	A
Chr4	TAIR8	9257157	substitution	C	G	***	Chr4	TAIR9	9257135		substitution	C	G
Chr4	TAIR8	9258137	deletion	C		***	Chr4	TAIR9	9258115	9258116	deletion	C
Chr4	TAIR8	9259264	substitution	C	G	***	Chr4	TAIR9	9259241		substitution	C	G
Chr4	TAIR8	9260118	substitution	C	T	***	Chr4	TAIR9	9260095		substitution	C	T
Chr4	TAIR8	9271022	insertion	C		***	Chr4	TAIR9	9270999		insertion	C
Chr4	TAIR8	9280168	substitution	A	C	***	Chr4	TAIR9	9280146		substitution	A	C
Chr4	TAIR8	9291817	substitution	A	G	***	Chr4	TAIR9	9291795		substitution	A	G
Chr4	TAIR8	9291823	substitution	A	G	***	Chr4	TAIR9	9291801		substitution	A	G
Chr4	TAIR8	9294797	deletion	G		***	Chr4	TAIR9	9294775	9294776	deletion	G
Chr4	TAIR8	9298235	substitution	T	A	***	Chr4	TAIR9	9298212		substitution	T	A
Chr4	TAIR8	9313829	substitution	C	A	***	Chr4	TAIR9	9313806		substitution	C	A
Chr4	TAIR8	9314183	insertion	C		***	Chr4	TAIR9	9314160		insertion	C
Chr4	TAIR8	9317432	substitution	C	G	***	Chr4	TAIR9	9317410		substitution	C	G
Chr4	TAIR8	9329709	substitution	G	A	***	Chr4	TAIR9	9329687		substitution	G	A
Chr4	TAIR8	9425005	insertion	G		***	Chr4	TAIR9	9424983		insertion	G
Chr4	TAIR8	9431920	insertion	A		***	Chr4	TAIR9	9431899		insertion	A
Chr4	TAIR8	9434557	insertion	A		***	Chr4	TAIR9	9434537		insertion	A
Chr4	TAIR8	9438362	substitution	T	A	***	Chr4	TAIR9	9438343		substitution	T	A
Chr4	TAIR8	9594883	substitution	A	G	***	Chr4	TAIR9	9594864		substitution	A	G
Chr4	TAIR8	9598063	insertion	C		***	Chr4	TAIR9	9598044		insertion	C
Chr4	TAIR8	9598078	insertion	G		***	Chr4	TAIR9	9598060		insertion	G
Chr4	TAIR8	9601982	insertion	C		***	Chr4	TAIR9	9601965		insertion	C
Chr4	TAIR8	9608779	substitution	C	T	***	Chr4	TAIR9	9608763		substitution	C	T
Chr4	TAIR8	9610209	deletion	G		***	Chr4	TAIR9	9610193	9610194	deletion	G
Chr4	TAIR8	9610659	insertion	T		***	Chr4	TAIR9	9610642		insertion	T
Chr4	TAIR8	9610729	insertion	C		***	Chr4	TAIR9	9610713		insertion	C
Chr4	TAIR8	9611133	insertion	C		***	Chr4	TAIR9	9611118		insertion	C
Chr4	TAIR8	9611303	insertion	C		***	Chr4	TAIR9	9611289		insertion	C
Chr4	TAIR8	9611704	deletion	G		***	Chr4	TAIR9	9611691	9611692	deletion	G
Chr4	TAIR8	9612232	insertion	C		***	Chr4	TAIR9	9612218		insertion	C
Chr4	TAIR8	9612460	insertion	C		***	Chr4	TAIR9	9612447		insertion	C
Chr4	TAIR8	9614126	substitution	C	A	***	Chr4	TAIR9	9614114		substitution	C	A
Chr4	TAIR8	9614167	substitution	C	A	***	Chr4	TAIR9	9614155		substitution	C	A
Chr4	TAIR8	9614245	substitution	C	A	***	Chr4	TAIR9	9614233		substitution	C	A
Chr4	TAIR8	9615720	insertion	C		***	Chr4	TAIR9	9615708		insertion	C
Chr4	TAIR8	9661422	insertion	A		***	Chr4	TAIR9	9661411		insertion	A
Chr4	TAIR8	9661961	substitution	T	C	***	Chr4	TAIR9	9661951		substitution	T	C
Chr4	TAIR8	9662601	insertion	A		***	Chr4	TAIR9	9662591		insertion	A
Chr4	TAIR8	9666705	deletion	A		***	Chr4	TAIR9	9666696	9666697	deletion	A
Chr4	TAIR8	9666732	deletion	A		***	Chr4	TAIR9	9666722	9666723	deletion	A
Chr4	TAIR8	9667889	deletion	T		***	Chr4	TAIR9	9667878	9667879	deletion	T
Chr4	TAIR8	9667896	deletion	T		***	Chr4	TAIR9	9667884	9667885	deletion	T
Chr4	TAIR8	9668173	insertion	A		***	Chr4	TAIR9	9668160		insertion	A
Chr4	TAIR8	9672325	deletion	G		***	Chr4	TAIR9	9672313	9672314	deletion	G
Chr4	TAIR8	9672457	insertion	A		***	Chr4	TAIR9	9672444		insertion	A
Chr4	TAIR8	9672997	deletion	T		***	Chr4	TAIR9	9672985	9672986	deletion	T
Chr4	TAIR8	9673008	deletion	T		***	Chr4	TAIR9	9672995	9672996	deletion	T
Chr4	TAIR8	9673332	insertion	A		***	Chr4	TAIR9	9673318		insertion	A
Chr4	TAIR8	9679425	insertion	A		***	Chr4	TAIR9	9679412		insertion	A
Chr4	TAIR8	9680514	substitution	C	A	***	Chr4	TAIR9	9680502		substitution	C	A
Chr4	TAIR8	9680635	substitution	C	A	***	Chr4	TAIR9	9680623		substitution	C	A
Chr4	TAIR8	9722677	substitution	C	G	***	Chr4	TAIR9	9722665		substitution	C	G
Chr4	TAIR8	9759486	deletion	G		***	Chr4	TAIR9	9759474	9759475	deletion	G
Chr4	TAIR8	9951766	insertion	G		***	Chr4	TAIR9	9951753		insertion	G
Chr4	TAIR8	10018844	substitution	A	T	***	Chr4	TAIR9	10018832		substitution	A	T
Chr4	TAIR8	10027232	deletion	G		***	Chr4	TAIR9	10027220	10027221	deletion	G
Chr4	TAIR8	10027351	insertion	T		***	Chr4	TAIR9	10027338		insertion	T
Chr4	TAIR8	10028601	substitution	G	A	***	Chr4	TAIR9	10028589		substitution	G	A
Chr4	TAIR8	10028606	substitution	C	A	***	Chr4	TAIR9	10028594		substitution	C	A
Chr4	TAIR8	10032493	substitution	A	C	***	Chr4	TAIR9	10032481		substitution	A	C
Chr4	TAIR8	10033995	substitution	C	T	***	Chr4	TAIR9	10033983		substitution	C	T
Chr4	TAIR8	10035570	substitution	C	T	***	Chr4	TAIR9	10035558		substitution	C	T
Chr4	TAIR8	10041886	insertion	G		***	Chr4	TAIR9	10041874		insertion	G
Chr4	TAIR8	10041904	insertion	G		***	Chr4	TAIR9	10041893		insertion	G
Chr4	TAIR8	10041918	deletion	G		***	Chr4	TAIR9	10041908	10041909	deletion	G
Chr4	TAIR8	10041921	deletion	T		***	Chr4	TAIR9	10041910	10041911	deletion	T
Chr4	TAIR8	10044468	substitution	T	C	***	Chr4	TAIR9	10044456		substitution	T	C
Chr4	TAIR8	10048081	insertion	A		***	Chr4	TAIR9	10048069		insertion	A
Chr4	TAIR8	10053788	substitution	T	A	***	Chr4	TAIR9	10053777		substitution	T	A
Chr4	TAIR8	10271494	substitution	C	T	***	Chr4	TAIR9	10271483		substitution	C	T
Chr4	TAIR8	10273587	insertion	G		***	Chr4	TAIR9	10273576		insertion	G
Chr4	TAIR8	10276350	deletion	C		***	Chr4	TAIR9	10276340	10276341	deletion	C
Chr4	TAIR8	10320065	substitution	C	G	***	Chr4	TAIR9	10320054		substitution	C	G
Chr4	TAIR8	10627037	insertion	G		***	Chr4	TAIR9	10627026		insertion	G
Chr4	TAIR8	10669050	substitution	A	G	***	Chr4	TAIR9	10669040		substitution	A	G
Chr4	TAIR8	10669053	substitution	C	A	***	Chr4	TAIR9	10669043		substitution	C	A
Chr4	TAIR8	10939347	substitution	T	C	***	Chr4	TAIR9	10939337		substitution	T	C
Chr4	TAIR8	10971887	substitution	T	C	***	Chr4	TAIR9	10971877		substitution	T	C
Chr4	TAIR8	10971904	substitution	A	T	***	Chr4	TAIR9	10971894		substitution	A	T
Chr4	TAIR8	10979236	insertion	T		***	Chr4	TAIR9	10979226		insertion	T
Chr4	TAIR8	10979429	deletion	A		***	Chr4	TAIR9	10979420	10979421	deletion	A
Chr4	TAIR8	10979435	deletion	A		***	Chr4	TAIR9	10979425	10979426	deletion	A
Chr4	TAIR8	11277590	insertion	G		***	Chr4	TAIR9	11277579		insertion	G
Chr4	TAIR8	11340005	substitution	G	T	***	Chr4	TAIR9	11339995		substitution	G	T
Chr4	TAIR8	11346839	substitution	A	T	***	Chr4	TAIR9	11346829		substitution	A	T
Chr4	TAIR8	11348715	deletion	T		***	Chr4	TAIR9	11348705	11348706	deletion	T
Chr4	TAIR8	11348801	substitution	G	T	***	Chr4	TAIR9	11348790		substitution	G	T
Chr4	TAIR8	11423677	substitution	T	C	***	Chr4	TAIR9	11423666		substitution	T	C
Chr4	TAIR8	11433344	deletion	G		***	Chr4	TAIR9	11433333	11433334	deletion	G
Chr4	TAIR8	11438304	substitution	T	C	***	Chr4	TAIR9	11438292		substitution	T	C
Chr4	TAIR8	11440335	substitution	C	A	***	Chr4	TAIR9	11440323		substitution	C	A
Chr4	TAIR8	11454794	substitution	C	T	***	Chr4	TAIR9	11454782		substitution	C	T
Chr4	TAIR8	11458337	substitution	G	A	***	Chr4	TAIR9	11458325		substitution	G	A
Chr4	TAIR8	11462753	substitution	T	C	***	Chr4	TAIR9	11462741		substitution	T	C
Chr4	TAIR8	11846366	substitution	G	A	***	Chr4	TAIR9	11846354		substitution	G	A
Chr4	TAIR8	11852534	insertion	G		***	Chr4	TAIR9	11852522		insertion	G
Chr4	TAIR8	12156576	insertion	C		***	Chr4	TAIR9	12156565		insertion	C
Chr4	TAIR8	13014968	substitution	G	T	***	Chr4	TAIR9	13014958		substitution	G	T
Chr4	TAIR8	13028886	insertion	C		***	Chr4	TAIR9	13028876		insertion	C
Chr4	TAIR8	13035522	substitution	T	C	***	Chr4	TAIR9	13035513		substitution	T	C
Chr4	TAIR8	13350915	insertion	C		***	Chr4	TAIR9	13350906		insertion	C
Chr4	TAIR8	13369945	substitution	C	G	***	Chr4	TAIR9	13369937		substitution	C	G
Chr4	TAIR8	13384787	insertion	G		***	Chr4	TAIR9	13384779		insertion	G
Chr4	TAIR8	13395245	deletion	A		***	Chr4	TAIR9	13395238	13395239	deletion	A
Chr4	TAIR8	13396553	insertion	G		***	Chr4	TAIR9	13396545		insertion	G
Chr4	TAIR8	13396624	insertion	G		***	Chr4	TAIR9	13396617		insertion	G
Chr4	TAIR8	13483949	substitution	T	A	***	Chr4	TAIR9	13483943		substitution	T	A
Chr4	TAIR8	13496990	substitution	A	T	***	Chr4	TAIR9	13496984		substitution	A	T
Chr4	TAIR8	13504606	substitution	C	G	***	Chr4	TAIR9	13504600		substitution	C	G
Chr4	TAIR8	13554554	substitution	A	G	***	Chr4	TAIR9	13554548		substitution	A	G
Chr4	TAIR8	13649671	substitution	G	T	***	Chr4	TAIR9	13649665		substitution	G	T
Chr4	TAIR8	13660159	substitution	G	A	***	Chr4	TAIR9	13660153		substitution	G	A
Chr4	TAIR8	14064437	deletion	C		***	Chr4	TAIR9	14064431	14064432	deletion	C
Chr4	TAIR8	14067725	insertion	G		***	Chr4	TAIR9	14067718		insertion	G
Chr4	TAIR8	15165547	substitution	A	T	***	Chr4	TAIR9	15165541		substitution	A	T
Chr4	TAIR8	15741462	insertion	G		***	Chr4	TAIR9	15741456		insertion	G
Chr4	TAIR8	15808872	insertion	T		***	Chr4	TAIR9	15808867		insertion	T
Chr4	TAIR8	15870964	substitution	C	A	***	Chr4	TAIR9	15870960		substitution	C	A
Chr4	TAIR8	15879560	insertion	G		***	Chr4	TAIR9	15879556		insertion	G
Chr4	TAIR8	15906626	substitution	A	G	***	Chr4	TAIR9	15906623		substitution	A	G
Chr4	TAIR8	15912361	substitution	G	T	***	Chr4	TAIR9	15912358		substitution	G	T
Chr4	TAIR8	15921748	substitution	A	G	***	Chr4	TAIR9	15921745		substitution	A	G
Chr4	TAIR8	15922634	substitution	A	G	***	Chr4	TAIR9	15922631		substitution	A	G
Chr4	TAIR8	15994043	substitution	T	A	***	Chr4	TAIR9	15994040		substitution	T	A
Chr4	TAIR8	16022561	deletion	C		***	Chr4	TAIR9	16022558	16022559	deletion	C
Chr4	TAIR8	16036255	substitution	G	T	***	Chr4	TAIR9	16036251		substitution	G	T
Chr4	TAIR8	16036272	substitution	G	T	***	Chr4	TAIR9	16036268		substitution	G	T
Chr4	TAIR8	16067636	insertion	C		***	Chr4	TAIR9	16067632		insertion	C
Chr4	TAIR8	16137173	insertion	A		***	Chr4	TAIR9	16137170		insertion	A
Chr4	TAIR8	16141689	insertion	T		***	Chr4	TAIR9	16141687		insertion	T
Chr4	TAIR8	16145065	insertion	A		***	Chr4	TAIR9	16145064		insertion	A
Chr4	TAIR8	16145427	deletion	G		***	Chr4	TAIR9	16145427	16145428	deletion	G
Chr4	TAIR8	16145436	deletion	AT		***	Chr4	TAIR9	16145435	16145436	deletion	AT
Chr4	TAIR8	16149592	substitution	T	A	***	Chr4	TAIR9	16149589		substitution	T	A
Chr4	TAIR8	16188039	insertion	C		***	Chr4	TAIR9	16188036		insertion	C
Chr4	TAIR8	16207825	deletion	T		***	Chr4	TAIR9	16207823	16207824	deletion	T
Chr4	TAIR8	16210148	deletion	T		***	Chr4	TAIR9	16210145	16210146	deletion	T
Chr4	TAIR8	16210222	deletion	T		***	Chr4	TAIR9	16210218	16210219	deletion	T
Chr4	TAIR8	16337356	insertion	G		***	Chr4	TAIR9	16337351		insertion	G
Chr4	TAIR8	16428227	substitution	G	A	***	Chr4	TAIR9	16428223		substitution	G	A
Chr4	TAIR8	16484721	deletion	A		***	Chr4	TAIR9	16484717	16484718	deletion	A
Chr4	TAIR8	16533917	substitution	T	G	***	Chr4	TAIR9	16533912		substitution	T	G
Chr4	TAIR8	16701091	substitution	C	A	***	Chr4	TAIR9	16701086		substitution	C	A
Chr4	TAIR8	16953990	insertion	G		***	Chr4	TAIR9	16953985		insertion	G
Chr4	TAIR8	16981615	insertion	G		***	Chr4	TAIR9	16981611		insertion	G
Chr4	TAIR8	17247695	deletion	C		***	Chr4	TAIR9	17247692	17247693	deletion	C
Chr4	TAIR8	17247781	deletion	A		***	Chr4	TAIR9	17247777	17247778	deletion	A
Chr4	TAIR8	17258726	substitution	A	T	***	Chr4	TAIR9	17258721		substitution	A	T
Chr4	TAIR8	17263864	deletion	G		***	Chr4	TAIR9	17263859	17263860	deletion	G
Chr4	TAIR8	17264859	substitution	T	C	***	Chr4	TAIR9	17264853		substitution	T	C
Chr4	TAIR8	17273615	insertion	T		***	Chr4	TAIR9	17273609		insertion	T
Chr4	TAIR8	17274261	insertion	G		***	Chr4	TAIR9	17274256		insertion	G
Chr4	TAIR8	17276585	deletion	A		***	Chr4	TAIR9	17276581	17276582	deletion	A
Chr4	TAIR8	17283711	substitution	G	C	***	Chr4	TAIR9	17283706		substitution	G	C
Chr4	TAIR8	17286716	insertion	A		***	Chr4	TAIR9	17286711		insertion	A
Chr4	TAIR8	17289242	substitution	G	A	***	Chr4	TAIR9	17289238		substitution	G	A
Chr4	TAIR8	17307382	substitution	A	T	***	Chr4	TAIR9	17307378		substitution	A	T
Chr4	TAIR8	17324609	deletion	T		***	Chr4	TAIR9	17324605	17324606	deletion	T
Chr4	TAIR8	17326697	deletion	C		***	Chr4	TAIR9	17326692	17326693	deletion	C
Chr4	TAIR8	17330774	deletion	C		***	Chr4	TAIR9	17330768	17330769	deletion	C
Chr4	TAIR8	17330793	insertion	A		***	Chr4	TAIR9	17330786		insertion	A
Chr4	TAIR8	17333539	insertion	T		***	Chr4	TAIR9	17333533		insertion	T
Chr4	TAIR8	17341734	insertion	A		***	Chr4	TAIR9	17341729		insertion	A
Chr4	TAIR8	17351527	insertion	A		***	Chr4	TAIR9	17351523		insertion	A
Chr4	TAIR8	17354494	insertion	G		***	Chr4	TAIR9	17354491		insertion	G
Chr4	TAIR8	17357683	substitution	C	A	***	Chr4	TAIR9	17357681		substitution	C	A
Chr4	TAIR8	17358578	insertion	C		***	Chr4	TAIR9	17358576		insertion	C
Chr4	TAIR8	17358583	insertion	A		***	Chr4	TAIR9	17358582		insertion	A
Chr4	TAIR8	17359511	insertion	T		***	Chr4	TAIR9	17359511		insertion	T
Chr4	TAIR8	17359552	deletion	A		***	Chr4	TAIR9	17359553	17359554	deletion	A
Chr4	TAIR8	17359825	substitution	A	C	***	Chr4	TAIR9	17359825		substitution	A	C
Chr4	TAIR8	17360258	insertion	A		***	Chr4	TAIR9	17360258		insertion	A
Chr4	TAIR8	17370694	substitution	G	A	***	Chr4	TAIR9	17370695		substitution	G	A
Chr4	TAIR8	17371743	insertion	TA		***	Chr4	TAIR9	17371744		insertion	TA
Chr4	TAIR8	17371744	insertion	G		***	Chr4	TAIR9	17371747		insertion	G
Chr4	TAIR8	17373398	deletion	T		***	Chr4	TAIR9	17373402	17373403	deletion	T
Chr4	TAIR8	17383209	substitution	T	A	***	Chr4	TAIR9	17383212		substitution	T	A
Chr4	TAIR8	17384770	substitution	A	C	***	Chr4	TAIR9	17384773		substitution	A	C
Chr4	TAIR8	17397852	substitution	C	G	***	Chr4	TAIR9	17397855		substitution	C	G
Chr4	TAIR8	17413659	substitution	T	C	***	Chr4	TAIR9	17413662		substitution	T	C
Chr4	TAIR8	17420206	insertion	C		***	Chr4	TAIR9	17420209		insertion	C
Chr4	TAIR8	17426250	substitution	A	C	***	Chr4	TAIR9	17426254		substitution	A	C
Chr4	TAIR8	17428027	insertion	G		***	Chr4	TAIR9	17428031		insertion	G
Chr4	TAIR8	17434023	insertion	G		***	Chr4	TAIR9	17434028		insertion	G
Chr4	TAIR8	18458057	substitution	T	A	***	Chr4	TAIR9	18458063		substitution	T	A
Chr4	TAIR8	18578130	substitution	A	G	***	Chr4	TAIR9	18578136		substitution	A	G
Chr4	TAIR8	18580947	insertion	AA		***	Chr4	TAIR9	18580953		insertion	AA
Chr4	TAIR8	18581145	insertion	CT		***	Chr4	TAIR9	18581153		insertion	CT
Chr4	TAIR8	18581342	insertion	GG		***	Chr4	TAIR9	18581352		insertion	GG
Chr4	TAIR8	18581539	insertion	AA		***	Chr4	TAIR9	18581551		insertion	AA
Chr5	TAIR8	3624	insertion	G		***	Chr5	TAIR9	3624		insertion	G
Chr5	TAIR8	225979	insertion	A		***	Chr5	TAIR9	225980		insertion	A
Chr5	TAIR8	225989	insertion	C		***	Chr5	TAIR9	225991		insertion	C
Chr5	TAIR8	426632	deletion	G		***	Chr5	TAIR9	426635	426636	deletion	G
Chr5	TAIR8	710905	deletion	T		***	Chr5	TAIR9	710907	710908	deletion	T
Chr5	TAIR8	1038587	substitution	T	G	***	Chr5	TAIR9	1038588		substitution	T	G
Chr5	TAIR8	1097774	deletion	G		***	Chr5	TAIR9	1097775	1097776	deletion	G
Chr5	TAIR8	1097867	deletion	C		***	Chr5	TAIR9	1097867	1097868	deletion	C
Chr5	TAIR8	2848847	substitution	A	C	***	Chr5	TAIR9	2848846		substitution	A	C
Chr5	TAIR8	4038336	deletion	T		***	Chr5	TAIR9	4038335	4038336	deletion	T
Chr5	TAIR8	4038480	deletion	C		***	Chr5	TAIR9	4038478	4038479	deletion	C
Chr5	TAIR8	4472484	insertion	G		***	Chr5	TAIR9	4472481		insertion	G
Chr5	TAIR8	4591928	deletion	A		***	Chr5	TAIR9	4591926	4591927	deletion	A
Chr5	TAIR8	4896803	substitution	A	G	***	Chr5	TAIR9	4896800		substitution	A	G
Chr5	TAIR8	5638930	insertion	C		***	Chr5	TAIR9	5638927		insertion	C
Chr5	TAIR8	7174421	deletion	A		***	Chr5	TAIR9	7174419	7174420	deletion	A
Chr5	TAIR8	8799791	substitution	A	C	***	Chr5	TAIR9	8799788		substitution	A	C
Chr5	TAIR8	9807548	deletion	T		***	Chr5	TAIR9	9807545	9807546	deletion	T
Chr5	TAIR8	11211679	substitution	G	A	***	Chr5	TAIR9	11211675		substitution	G	A
Chr5	TAIR8	11226526	deletion	Nx14108		***	Chr5	TAIR9	11226522	11226523	deletion	Nx14108
Chr5	TAIR8	11241928	substitution	G	C	***	Chr5	TAIR9	11227816		substitution	G	C
Chr5	TAIR8	11242113	substitution	C	G	***	Chr5	TAIR9	11228001		substitution	C	G
Chr5	TAIR8	11242399	deletion	Nx2134		***	Chr5	TAIR9	11228287	11228288	deletion	Nx2134
Chr5	TAIR8	11245005	deletion	Nx984		***	Chr5	TAIR9	11228759	11228760	deletion	Nx984
Chr5	TAIR8	11662967	substitution	C	T	***	Chr5	TAIR9	11645737		substitution	C	T
Chr5	TAIR8	11788160	substitution	A	G	***	Chr5	TAIR9	11770930		substitution	A	G
Chr5	TAIR8	11800919	substitution	A	C	***	Chr5	TAIR9	11783689		substitution	A	C
Chr5	TAIR8	11820771	substitution	A	G	***	Chr5	TAIR9	11803541		substitution	A	G
Chr5	TAIR8	15642817	insertion	A		***	Chr5	TAIR9	15625587		insertion	A
Chr5	TAIR8	15652856	insertion	T		***	Chr5	TAIR9	15635627		insertion	T
Chr5	TAIR8	17072077	substitution	C	A	***	Chr5	TAIR9	17054849		substitution	C	A
Chr5	TAIR8	17518010	insertion	T		***	Chr5	TAIR9	17500782		insertion	T
Chr5	TAIR8	18952108	substitution	C	T	***	Chr5	TAIR9	18934881		substitution	C	T
Chr5	TAIR8	19301885	insertion	C		***	Chr5	TAIR9	19284658		insertion	C
