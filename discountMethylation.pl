#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use 5.010_000;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use FastaReader;
use DBI;
use Log::Log4perl qw/:easy/;
use List::Util qw/max/;
use List::MoreUtils qw/all/;
use DBI;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use Module::Load::Conditional qw/can_load/;
use File::CountLines qw/count_lines/;
use File::Temp qw/mktemp tempfile/;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if !$opt_file || !$opt_output_prefix || !$opt_reference;
Log::Log4perl->easy_init({ level => $DEBUG, layout => '%d{HH:mm:ss} %.1p > %m%n' });
my $logger = get_logger();

#######################################################################
# progress bar!

my $has_progressbar = $opt_progress_bar && can_load(modules => { 'File::CountLines' => undef, 'Term::ProgressBar' => undef });
my $linecount;
my $pb_increment = 3000; # heuristically determined
my $pb;
my $last_pb_update=0;

if ($has_progressbar){
    $logger->info("Has Term::ProgressBar, using it");
    $linecount = count_lines($opt_file);
    $logger->info("$linecount lines total to process");
    $pb = Term::ProgressBar->new({count => $linecount,ETA => 'linear'});
    $pb->minor(1);
}

#######################################################################
# Database 

{
    my $counter = 0;
    my $increment = 10000;
    my $filename = $opt_memory ? ':memory:' : (tempfile("$opt_output_prefix.tmp.XXXXX", UNLINK => 1))[1];
    unlink $filename if -f $filename;
    my $dbh = DBI->connect("dbi:SQLite:dbname=$filename","","", {RaiseError => 1, AutoCommit => 0});
    $dbh->do("PRAGMA journal_mode = OFF");
    $dbh->do("PRAGMA cache_size = 80000");
    $dbh->do("create table methyl (seq, coord integer, context, c integer default 0, t integer default 0)");
    $dbh->do("create index idx on methyl (seq,coord)");

    my $checker = $dbh->prepare("select count(seq),context from methyl where seq=? and coord=?");

    my %updater = map {
        $_ => $dbh->prepare("update methyl set $_=(select $_ from methyl where seq=? and coord=?)+1 where seq=? and coord=?"),
    } qw/c t/;

	my %inserter = map {
        $_ => $dbh->prepare("insert into methyl(seq, coord, context, $_) values (?,?,?,1)")
    } qw/c t/;

    sub disconnect { $dbh->disconnect }

    # if record exists, return context
    sub record_exists{
        my ($seq, $coord) = @_;
        $checker->execute($seq,$coord);
        my $row = $checker->fetch();
        return $row->[0] ? $row->[1] : 0;
    }

    sub record_methylation{
        my ($seq, $coord, $context) = @_;
        $context = lc $context;
        my $base;

        # normalize tg/thg/thh => cg/chg/chh, but record the first base
        if ($context =~ s/^(c|t)/c/){
            $base = $1;
        }
        else {
            die "$context not c/t?";
        }

        #say STDERR join ",", @_;
        if (my $existing_context = record_exists($seq,$coord)){
            #say STDERR "record exists";
            die "incompatible context ($existing_context, $context) at same coord? bug!" if $existing_context ne $context;
            $updater{$base}->execute($seq,$coord,$seq,$coord);
        }
        else{
            #say STDERR "record d.n.exists";
            $inserter{$base}->execute($seq, $coord, $context);
        }
        if (++$counter % $increment == 0){
            $dbh->commit();
        }
    }

    sub record_output{
        # prefix-$seq.single-c-$context.gff.merged
        my $prefix = shift;

        $dbh->commit();

        my $sth = $dbh->prepare('select seq, coord, context, c, t from methyl order by seq, coord');
        $sth->execute();

        my %filehandles;
        my %temp2real;

        while (defined(my $row = $sth->fetch)){
            #say join "|", @$row;
            my ($seq, $coord, $context, $c, $t) = @$row;
            $context = uc $context;

            my $key = $context;
            if (! exists $filehandles{$key}){
                my $file = "$prefix.single-c.$context.gff.merged";
                my $tmpfile = mktemp($file . ".tmp.XXXX");
                open my $writer, q{>}, $tmpfile;
                $filehandles{$key} = $writer;
                $temp2real{$tmpfile} = $file;
            }
            die "why is \$c + \$t == 0? bug, dying" if $c + $t == 0;

            my $score = sprintf("%.4f", $c/($c+$t));
            say {$filehandles{$key}} join "\t", $seq, q{.}, $context, $coord, $coord, $score, q{.}, q{.}, "c=$c;t=$t";
        }

        close $_ for values %filehandles;

        while (my ($tmp,$real) = each %temp2real) {
            rename $tmp, $real;
        }
    }
}

#######################################################################

{
    my %stats;

    sub count_methylation{
        my ($gff_line, $fr) = @_;
        #die "$gff_line\n" . Dumper $sequence_lengths;

        my $filtered = $gff_line =~ s/\*$//;

        my @split = split /\t/, $gff_line;

        die "gff should be 9 columed..." if @split != 9;

        my ($seq, $start,$end,$strand) = @split[0,3,4,6];

        # sequence - create entry in stats if necessary.

        return if ($seq eq '.');
        if (! exists $stats{$seq}){
            if ($opt_dinucleotide){
                $stats{$seq} = {
                    bp         => 0,
                    unfiltered => { map {$_ => 0} qw/c cg cc ca ct t tg tc ta tt/ },
                    filtered   => { map {$_ => 0} qw/c cg cc ca ct t tg tc ta tt/ },
                };
            }
            else{
                $stats{$seq} = {
                    bp         => 0,
                    unfiltered => { map {$_ => 0} qw/c cg chh chg t tg thh thg/ },
                    filtered   => { map {$_ => 0} qw/c cg chh chg t tg thh thg/ },
                };
            }
        }

        my $unfiltered_count = $stats{$seq}{unfiltered}; # deref them here, once, for speed
        my $filtered_count   = $stats{$seq}{filtered};

        # Grab sequences

        my ($read_seq, $target_seq);
        #if ($split[2] =~ /([ATCGN]+)$/){
        if ($split[2] =~ /([A-Z]+)$/){
            $read_seq= $1;
        }
        #if ($split[8] =~ /target=([ATCGN]+)$/){
        if ($split[8] =~ /target=([A-Z]+)$/){
            $target_seq = $1;
        }
        if (!  defined $read_seq || ! defined $target_seq){
            $logger->debug("couldn't parse $gff_line ?");
            return;
        }
        die "can't find read or target seq"  unless (defined $read_seq && defined $target_seq);

        my @target_bases = split //,$target_seq;
        my @read_bases = (q{.}, q{.}, split(//, $read_seq), q{.}, q{.});

        # check length

        if (length($read_seq) != ($end-$start+1) || length $target_seq != 4 + length $read_seq){
            die "read size mismatch\n$read_seq\n$target_seq";
        }
        else{
            $stats{$seq}{bp} += length($read_seq);
        }

        # reverse?

        my $reverse = $strand eq '+' ? 0 : $strand eq '-' ? 1 : die "strand should be + or -";

        READ:
        for (my $strand_coord = $start; $strand_coord <= $end; ++$strand_coord){
            my $i = $strand_coord - $start + 2;
            my $abs_coord = $reverse ? $fr->get_length($split[0]) - $strand_coord + 1 : $strand_coord;
            my $context;

            # first position
            my $methylation;
            if ($target_bases[$i] eq 'C'){
                given ($read_bases[$i]){
                    when ('C'){ $methylation = 1; }
                    when ('T'){ $methylation = 0; }
                    default { next READ; }
                }
            }
            else{
                next READ;
            }

            # second/third position

            if ($opt_dinucleotide){
                given ($target_bases[$i+1]){
                    when ('G'){ $context = $methylation ? 'cg' : 'tg'; }
                    when ('A'){ $context = $methylation ? 'ca' : 'ta'; }
                    when ('C'){ $context = $methylation ? 'cc' : 'tc'; }
                    when ('T'){ $context = $methylation ? 'ct' : 'tt'; }
                    default { next READ; }
                }
            }
            else{
                if ($target_bases[$i + 1] eq 'G'){
                    $context = $methylation ? 'cg' : 'tg';
                }
                elsif ($target_bases[$i + 2] eq 'G'){
                    $context = $methylation ? 'chg' : 'thg';
                }
                else {
                    $context = $methylation ? 'chh' : 'thh';
                }
            }

            if ($filtered){
                ++$filtered_count->{$context};
                ++$filtered_count->{$methylation ? 'c' : 't'};
            } else{
                ++$unfiltered_count->{$context};
                ++$filtered_count->{$methylation ? 'c' : 't'};
            }

            record_methylation($seq,$abs_coord,$context);
        }
    }

    sub rat{
        my ($x,$y) = @_;
        if ($x+$y != 0){
            return sprintf("%.6f", $x/($x+$y));
        }
        else {
            return "0.000000";
        }
    }

    # this thing looks ridiculous...
    sub print_freq{
        my $prefix = shift;

        open my $out, '>', "$prefix.freq";

        my @output;

        if ($opt_dinucleotide){
            push @output, [qw/seq bp overlaps 
            C CG CT CA CC 
            T TG TT TA TC
            C_ratio CG_ratio CT_ratio CA_ratio CC_ratio 

            filtered_C filtered_CG filtered_CT filtered_CA filtered_CC 
            filtered_T filtered_TG filtered_TT filtered_TA filtered_TC

            filtered_C_ratio filtered_CG_ratio filtered_CT_ratio filtered_CA_ratio filtered_CC_ratio 
            /];
        }
        else{
            push @output, [qw/seq bp overlaps 
            C CG CHG CHH 
            T TG THG THH 
            C_ratio CG_ratio CHG_ratio CHH_ratio 
            filtered_C filtered_CG filtered_CHG filtered_CHH 
            filtered_T filtered_TG filtered_THG filtered_THH 
            filtered_C_ratio filtered_CG_ratio filtered_CHG_ratio filtered_CHH_ratio 
            /];
        }
        
        for my $seq (sort keys %stats){
            if ($opt_dinucleotide){
                push @output, [
                $seq, $stats{$seq}{bp}, 0,

                $stats{$seq}{unfiltered}{c}, $stats{$seq}{unfiltered}{cg}, $stats{$seq}{unfiltered}{ct},$stats{$seq}{unfiltered}{ca}, $stats{$seq}{unfiltered}{cc}, 
                $stats{$seq}{unfiltered}{t}, $stats{$seq}{unfiltered}{tg}, $stats{$seq}{unfiltered}{tt},$stats{$seq}{unfiltered}{ta}, $stats{$seq}{unfiltered}{tc}, 

                rat($stats{$seq}{unfiltered}{c}  ,$stats{$seq}{unfiltered}{t} ),
                rat($stats{$seq}{unfiltered}{cg} ,$stats{$seq}{unfiltered}{tg} ),
                rat($stats{$seq}{unfiltered}{ct} ,$stats{$seq}{unfiltered}{tt} ),
                rat($stats{$seq}{unfiltered}{ca} ,$stats{$seq}{unfiltered}{ta} ),
                rat($stats{$seq}{unfiltered}{cc} ,$stats{$seq}{unfiltered}{tc} ),

                $stats{$seq}{filtered}{c}, $stats{$seq}{filtered}{cg}, $stats{$seq}{filtered}{ct},$stats{$seq}{filtered}{ca}, $stats{$seq}{filtered}{cc}, 
                $stats{$seq}{filtered}{t}, $stats{$seq}{filtered}{tg}, $stats{$seq}{filtered}{tt},$stats{$seq}{filtered}{ta}, $stats{$seq}{filtered}{tc}, 

                rat($stats{$seq}{filtered}{c}  ,$stats{$seq}{filtered}{t} ),
                rat($stats{$seq}{filtered}{cg} ,$stats{$seq}{filtered}{tg} ),
                rat($stats{$seq}{filtered}{ct} ,$stats{$seq}{filtered}{tt} ),
                rat($stats{$seq}{filtered}{ca} ,$stats{$seq}{filtered}{ta} ),
                rat($stats{$seq}{filtered}{cc} ,$stats{$seq}{filtered}{tc} ),
                ]
            }
            else{
                push @output, [
                $seq, $stats{$seq}{bp}, 0,

                $stats{$seq}{unfiltered}{c}, $stats{$seq}{unfiltered}{cg}, $stats{$seq}{unfiltered}{chg}, $stats{$seq}{unfiltered}{chh},
                $stats{$seq}{unfiltered}{t}, $stats{$seq}{unfiltered}{tg}, $stats{$seq}{unfiltered}{thg}, $stats{$seq}{unfiltered}{thh},

                rat($stats{$seq}{unfiltered}{c}  ,$stats{$seq}{unfiltered}{t} ),
                rat($stats{$seq}{unfiltered}{cg} ,$stats{$seq}{unfiltered}{tg} ),
                rat($stats{$seq}{unfiltered}{chg},$stats{$seq}{unfiltered}{thg} ),
                rat($stats{$seq}{unfiltered}{chh},$stats{$seq}{unfiltered}{thh} ),

                $stats{$seq}{filtered}{c}, $stats{$seq}{filtered}{cg}, $stats{$seq}{filtered}{chg}, $stats{$seq}{filtered}{chh},
                $stats{$seq}{filtered}{t}, $stats{$seq}{filtered}{tg}, $stats{$seq}{filtered}{thg}, $stats{$seq}{filtered}{thh},

                rat($stats{$seq}{filtered}{c}  , $stats{$seq}{filtered}{t} ),
                rat($stats{$seq}{filtered}{cg} , $stats{$seq}{filtered}{tg} ),
                rat($stats{$seq}{filtered}{chg}, $stats{$seq}{filtered}{thg} ),
                rat($stats{$seq}{filtered}{chh}, $stats{$seq}{filtered}{thh} ),
                ]
            }
        }

        #say scalar(@$_) for @output;
        my $numcols = $opt_dinucleotide ? 33 : 27;

        die "uneven number of lines in freq? dying" unless all { $numcols == scalar @$_ } @output;

        for my $i (0..$numcols-1) {
            my $line = join "\t", map { $_->[$i] } @output;
            $logger->debug($line);
            say $out $line;
        }
        close $out;
    }

}

#######################################################################
# Main body

my $fr = FastaReader->new(file => $opt_reference, normalize => 0);

open my $in, '<', $opt_file;

while (defined(my $line = <$in>)){
    #chomp $line;
    $line =~ tr/\n\r//d;
    count_methylation($line, $fr);

    if ($. % $pb_increment == 0){
        if ($has_progressbar){
            $pb->update($.);
        }
        else{
            $logger->debug($.);
        }
    }
}
close $in;


$pb->update($linecount) if $has_progressbar;
$logger->info("Done processing! creating single-c and freq file");

record_output($opt_output_prefix);
print_freq($opt_output_prefix);
disconnect;

#######################################################################
# DONE

=head1 NAME

discountMethylation.pl - ...

=head1 SYNOPSIS

Usage examples:

 discountMethylation.pl [options]...

=head1 OPTIONS

=over

=item  -o <prefix> | --output-prefix <prefix>

=for Euclid
    prefix.default:     '-'

=item  <file> 

=for Euclid
    file.type:        readable

=item  -r <fasta> | --reference <fasta>

=for Euclid
    fasta.type:        readable

=item  -pb | --progress-bar 

=item  -d | --dinucleotide 

=item --help | -h

=item  -m | --memory 

=back

=cut
__DATA__
chr3	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:1804:1117#0/1:TTTTATGATGTGGTAATTTATTATTGGATGGGAAGTTTGAT	1567476	1567516	1	+	0	target=TGTCCCATGACGTGGCAACCTATTACTGGATGGGAAGTTCGACCG
chr1	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:1996:1117#0/1:TGGATTGTAAGAAATAGAGAAGAAGTTAATTAGTGATTGAA	17479814	17479854	1	-	0	target=TGTGGACTGCAAGAAATAGAGAAGAAGCTAATTAGTGACTGAATA
chr5	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:2211:1114#0/1:AGTTTATTATTGAGAGAGGTTTGTGGAAATATGAGAGAATG	95815	95855	1	-	0	target=TTAGTTTATTACTGAGAGAGGCTTGTGGAAACACGAGAGAATGCG
chrc	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:2375:1110#0/1:TTTGAATAGAGGTTATGGATTTTNTTTTTGTAGAAGTAATT	53181	53221	1	+	1	target=TGCTTGAATAGAGGTTATGGACCCCTTTTTCGTAGAAGTAATTCT
.	NM/NM	SOLEXA2_0531_FC62W31AAXX:4:1:2392:1116#0/1:TTTTGGTTATATTATGAAAGTTTTATGAAGTTATAAGAAGT	0	0	0	.	.	.
chr5	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:3132:1115#0/1:TGTGTAAGATATTTTTTGAAAGAGAATTTTTTTAGAGATGA	14708724	14708764	1	+	0	target=AGCGCGTAAGATATTTCTTGAAAGAGAATTTTTCTAGAGATGATT
chrc	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:3177:1117#0/1:TTTATTTGATATGGGTAGGAGAAATGGTATTTTTTTTTTAA	138223	138263	1	+	0	target=GTTCTATTCGATACGGGTAGGAGAAACGGTATTCTTTTCTTAAAC
chr3	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:3197:1115#0/1:AAATAAGATATTTTTATTTATGTTTGTTGTGTGTAAGAAAA	6073818	6073858	1	-	0	target=TTAAACAAGACATCTCTATCCACGCTTGCCGTGCGTAAGAAAAGG
chrm	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:3227:1111#0/1:TGTTTTGGTAGGAAGTATTTTTTNTGAAGAATTGAATAAAA	148388	148428	1	+	1	target=CTCGCTTTGGCAGGAAGTATTTTTTTTGAAGAATTGAACAAAAAA
chrm	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:3470:1115#0/1:AATTTAAGATTTTTTAGAGTATTTTTGTTTTATTAATTTGT	42318	42358	1	-	0	target=ACAATCCAAGACCCCTCAGAGTATTTCTGTTTTATCAATTCGTTT
chr5	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:3700:1110#0/1:GTAAGAATAGTTTATATGAATGGNTTTTTTATTTTAAGATT	3167015	3167055	1	-	1	target=TTGCAAGAACAGTTTATACGAACGGTCTCTTTATTTCAAGATTTT
chr4	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:4405:1109#0/1:GAATTCGTGTATATGTATTAATTNTCGGGATTTCGTGATTT	1991034	1991074	1	+	1	target=AAGAACCCGTGTATATGCATCAACCTTCGGGATTCCGTGATCCCA
chr3	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:4434:1110#0/1:ATATGTATTAGTTGTAAAATATGNGGTTAATTTTATATTAG	7191215	7191255	1	-	1	target=TAATATGTATCAGCTGCAAAATATGAGGTTAATTTTACATCAGTT
chr1	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:4556:1111#0/1:ATTATGATTTTTGATGTGGTATANAATAATATAATAAATTA	21622283	21622323	1	+	1	target=AAATTATGATCCTTGATGTGGTATAGAACAACATAATAAATTAAT
chr5	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:4638:1110#0/1:AATTTGATAATTTATTAAGTTTTNAATTTAAATAAGTAAAT	8436366	8436406	1	+	1	target=GTAATTCGACAATTCATCAAGTCTCAAATCTAAATAAGTAAACAT
chr2	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:4681:1116#0/1:AGAAGAAAGGTGAAGAAGTTTTTAGTTTTGGTTTTGATTAA	11634907	11634947	1	+	0	target=TGAGAAGAAAGGTGAAGAAGCTCCTAGTCTTGGTCTTGATCAAGT
chr3	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:4813:1110#0/1:ATTTAATATGTGATGAAAGTATGNAAGTTTTGTTTTTTTGG	230904	230944	1	+	1	target=TAACTTAATATGCGACGAAAGCATGGAAGCTTCGTCTCTTTGGGA
chrc	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:5070:1117#0/1:TATTATTGTGGTAAAGGTTGTAATGTTAGAGGAATAATTAT	152908	152948	1	+	0	target=AGCATCATTGTGGTAAAGGTCGTAATGCCAGAGGAATAATTACCG
chrc	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:5111:1110#0/1:AGTATTTTTTGTTGTGGTTTGGGNAGTAAAGGGTGTTTTTT	79140	79180	1	+	1	target=ATAGCATTTCCTGCTGCGGTTTGGGCAGCAAAGGGTGTTCCTCTT
chr2	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:5164:1112#0/1:GAATTTTTAATTAGAAGATGTTTTGTATGTTTTGATGTAAT	12407634	12407674	1	+	0	target=GAGAATCTTCAATCAGAAGATGTCCCGCATGTTTTGATGTAACCT
chr5	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:5203:1113#0/1:TTTGTTATTATATTGTTTAGGTATTTTTTGAATAAAATTAT	16466392	16466432	1	+	0	target=AGCCCGCCACCACACTGCTCAGGTACCTTCTGAATAAAATTACTA
.	NM/NM	SOLEXA2_0531_FC62W31AAXX:4:1:5388:1114#0/1:TTTAAAATATTAATTAATTTTTTTTTTGTTTTTTAAAGTTT	0	0	0	.	.	.
.	NM/NM	SOLEXA2_0531_FC62W31AAXX:4:1:5444:1109#0/1:AAGATATAAAGTTAAAGATTTATNTGGATTTTGGTTATATT	0	0	0	.	.	.
chr4	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:5492:1112#0/1:TTAAAGTGGTTGTTGTGGTTGTTAGGTTTTAAGGTTAAGGT	14379089	14379129	1	+	0	target=GGCTAAAGCGGCTGCTGCGGTTGTTAGGCTCCAAGGTCAAGGCAA
chrc	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:5644:1109#0/1:GATTGTAAATTATAAGATTTTTTNAATTGTTTTTTTTGAAA	3176	3216	1	+	1	target=TTGATTGTAAATTATAAGATTTTTTTAATTGTTTTCCTTGAAAAG
chr5	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:5744:1114#0/1:AAACGAAGATTATATTTATTTGGAGTTGTATGGATTATCTG	13175521	13175561	1	+	0	target=ATAAACGAAGATTATATTCACTTGGAGCCGTATGGATTACCTGAA
chr1	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:5987:1114#0/1:TTTATTTATTTTCGATTATTGATATTTTTTAATTGTGTTAT	26289299	26289339	1	+	0	target=CGTTTATCTATTTTCGATTACTGATACTTTCTAACCGTGTTATGT
.	NM/NM	SOLEXA2_0531_FC62W31AAXX:4:1:6130:1111#0/1:TGGGTGATGAAGTAGGTTGTTGTNTATTTGATTTTGTTTAT	0	0	0	.	.	.
chr5	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:6323:1112#0/1:ATTATTTTTTTGGTTTGGTGATATAATATTATTATTTATGA	16485535	16485575	1	+	0	target=TAATTATTTCTTCGGTTTGGTGATACAACATTACTATTTATGATT
chr3	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:6363:1114#0/1:TTATTGATATTTTTATTATTTTTATTTTTAAAATGTAAGTG	15611785	15611825	1	+	0	target=ACCTATCGATATCTCTATTATTCTCATCTTCAAAACGTAAGTGTT
chr3	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:6449:1114#0/1:ACGATTTTTATTTTTGGTCGATAAATTTTCGGATTCGGTAA	9255571	9255611	1	-	0	target=CGACGATTTCCACTCCTGGTCGACAAATCCTCGGACCCGGTAAAC
.	NM/NM	SOLEXA2_0531_FC62W31AAXX:4:1:6591:1110#0/1:TTGATGGTGTAGTCAAAGTCCATNTGAGTTTTTATTTTTGT	0	0	0	.	.	.
chrc	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:6718:1116#0/1:AAAGATTTTTTTTTAAATTGTATATATGATTTTTATTGAAT	13879	13919	1	-	0	target=TGAAAGATCTCCCTCCAAACCGTACATACGACTTTCATCGAATAC
chr2	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:6978:1114#0/1:TTTTTTTTGAGTTTGGGTTATGGATGAAGGAGAAGAAGGAA	3187168	3187208	1	+	0	target=TTCCTTCTCCGAGCTCGGGTTATGGACGAAGGAGAAGAAGGAACC
chr2	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:7290:1114#0/1:TTTGTTTGGTTTGGAATATGGTTTAGAGAAATATTTAATTT	17180494	17180534	1	+	0	target=TCTTTGTTTGGTTTGGAACATGGTCTAGAGAAACATTCAATCCAG
.	NM/NM	SOLEXA2_0531_FC62W31AAXX:4:1:7436:1110#0/1:AGAGATTAGTTTATTTTTTGATTNGAAAAATAAGTTAATTT	0	0	0	.	.	.
.	NM/NM	SOLEXA2_0531_FC62W31AAXX:4:1:7501:1117#0/1:GAAGAGCGGTTCAGCAGGAATGCCGAGAACCGATCTCGTAT	0	0	0	.	.	.
.	NM/NM	SOLEXA2_0531_FC62W31AAXX:4:1:7535:1111#0/1:TAGGGGGGTTATTAGTTTATTGANTTTATTAAGTGTTATTT	0	0	0	.	.	.
chrc	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:7552:1117#0/1:AGGAATGGAAAAAATTGTAGAAAATTGAGTAATTATATAAT	31433	31473	1	-	0	target=AGAGGAATGGAAAAAATTGCAGAAAACCGAGCAATTATACAATAT
chr5	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:7636:1116#0/1:TATGAGTTATATTATGTAAATGTTAAATGTTTTTGGTGAGT	20144461	20144501	1	-	0	target=GCTATGAGCCACATCATGCAAATGTCAAATGTTCTCGGCGAGCAA
chrc	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:7679:1111#0/1:TATGAATAATTTTTTTTGGTATTNTATGTATATTTTTATGT	79283	79323	1	+	1	target=AACATGAATAACTCCCTTTGGTATTCTACGTACATTTTTACGTGA
chr3	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:7785:1111#0/1:TTAGTAATTAATTATATTTAAAANTTGAGTTAATATAATGT	10222308	10222348	1	+	1	target=GATTAGTAATTAATTACACTTAAAAATTGAGTTAACATAATGTAG
chrc	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:8104:1114#0/1:AGAATTTAGAATGTAAGAATGAATGGGTTTGTTTTGATTGT	45488	45528	1	-	0	target=TAAGAATTCAGAATGTAAGAATGAACGGGTTCGCTTTGACCGTTA
chr1	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:8155:1110#0/1:AGAAGGGATTTATAAATTAAATTNAATCGTCTAGTATACTG	15515317	15515357	1	+	1	target=CCAGAAGGGATCTATAAATCAAACTCAACCGCCTAGCATACTGAG
chr2	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:8268:1117#0/1:ATATTAGTAATATGGTTTAAAAATGATTTTGGTTTATATTA	15629450	15629490	1	-	0	target=CAACATCAGCAATATGGTTCAAAAATGATCCCGGTCTACATCAGT
chrm	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:8437:1111#0/1:ATTTGTATGTGTTAAGTATATAGNTAGTGTTTTTTTTGAGT	363215	363255	1	+	1	target=CGACTTGCATGTGTTAAGCATATAGCTAGCGTTCCTTCTGAGCCA
chrc	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:8513:1110#0/1:AAAGGAGTTGAATGAAATTAAAGNTTTATGTTTGGTTTTGA	9294	9334	1	+	1	target=ATAAAGGAGCCGAATGAAACCAAAGTTTCATGTTCGGTTTTGAAT
chr1	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:8543:1114#0/1:GAAAAAATAATTTGATTTATAAGGTAAATTTTAAGGTTTGA	436602	436642	1	+	0	target=GGGAAAAAATAACTCGACTTACAAGGTAAATCCCAAGGCTCGAAG
chrm	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:8677:1117#0/1:TGAGTTTTGAGTAGGTAGAGTTAATGTAAGTAGGGTTAATG	299692	299732	1	-	0	target=AGCGAGCTCTGAGCAGGCAGAGCTAATGCAAGTAGGGCCAACGAG
chr3	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:8781:1111#0/1:ATAATTTTTATCGTATTTATATANTTTTTTGCGTATTTAAT	6077478	6077518	1	-	1	target=GTATAATCTCTATCGTATTCACATAACTTTTTGCGTATTTAATCC
chr3	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:9077:1118#0/1:AATTTGTAAGAAAAATTAATTTGTTAATTTTTTTTGAATAT	13437933	13437973	1	+	0	target=CAAACTCGTAAGAAAAATTAACTCGTTAATTTTTTTTGAATATTT
chrc	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:9177:1112#0/1:GAGAAGGAAAAGGAGTTTTTTTGGTGGATATTTTTTGGTTT	128416	128456	1	+	0	target=AAGAGAAGGAAAAGGAGCTTCTTCGGTGGATACCTCTTGGTCCTG
chr5	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:9371:1117#0/1:GAAGAAATAAAGTTAAAGATTTATATGGATTTTGGTTATAT	11743743	11743783	1	+	0	target=TAGAAGAAACAAAGCCAAAGACTCATATGGACTTTGGCTACACCA
chrc	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:9445:1113#0/1:TTTTAGAATTATATGAAGTTAGTAAGTAAATAGTGAATTTA	130857	130897	1	-	0	target=GTTTTCAGAATTATATGAAGCCAGTAAGCAAACAGCGAATCCATG
.	NM/NM	SOLEXA2_0531_FC62W31AAXX:4:1:9732:1111#0/1:TGTTTAGTATATTGGTGATTTTTNTTATATGAAATGGGAGT	0	0	0	.	.	.
.	NM/NM	SOLEXA2_0531_FC62W31AAXX:4:1:10053:1113#0/1:TATGGATTTTGGTTATATTATGAAGGTTTTGAGAAGCTAGA	0	0	0	.	.	.
chrc	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:10217:1114#0/1:GTTAAAAATATAAAAATGATTAATTGATGTTGATTAAAATT	9344	9384	1	+	0	target=ACGTTAAAAATATAAAAATGATCAATCGACGTCGACTAAAACCCT
chr1	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:10457:1112#0/1:AGTTTTGTTTTTTTTGTTAATTCGAGTAATTTTATTTGTAT	27134077	27134117	1	-	0	target=TAAGCTTTGTCTTCCTTGCTAACCCGAGTAATCTTATTTGTACAT
chr5	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:10704:1112#0/1:AAATTTTACGGGTTATGTTTCGAGAAGAAAATAAAGAGAGA	10812243	10812283	1	+	0	target=CAAAATCTTACGGGCTATGCTTCGAGAAGAAAACAAAGAGAGATC
chr4	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:10810:1114#0/1:AGAATTATAGAATTGTTGTTGTGTTGATTGTTTTTAATGGA	7835324	7835364	1	-	0	target=GTAGAACTACAGAATTGTCGTCGTGTCGATCGTCTTCAATGGATT
chrc	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:11475:1115#0/1:TAGATAATTGAGGTTAATTTAGTTTTTAGATGAGTTTGGAT	102146	102186	1	-	0	target=GACAGACAATTGAGGCTAATCTAGCTCTCAGACGAGCTCGGACAC
chr3	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:11581:1114#0/1:TGTTAATGTTTTTTTTAGTATTTGGTTGAGATTTGTTTAAG	2174475	2174515	1	+	0	target=CGTGCCAATGTCTTCTCTAGTATCTGGCTGAGACTCGTTCAAGAA
chr2	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:11646:1113#0/1:TGGTAGGTTAGTTGGGTTATGAGTGTGGTGGGAATGGGTTT	16262319	16262359	1	-	0	target=GCTGGCAGGCCAGCCGGGCCATGAGCGCGGTGGGAACGGGCTTCC
chr5	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:11745:1111#0/1:AAGGGAATGGATGAGGTTTACAANAAAAGAATGGATTGGTA	14615223	14615263	1	-	1	target=GAAAGGGAATGGATGAGGCTTACAAGAAAAGAATGGATTGGTATC
chrc	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:11948:1116#0/1:TTAGAGATTATGAGTGTAATAGGAGTATTTGTTGATAAAAG	58392	58432	1	-	0	target=CTTCAGAGACTACGAGTGTAATAGGAGCATCCGTCGACAAAAGGA
chrc	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:12200:1117#0/1:AAGAATATTTTTTTTATTGGGATGATAATTTTGAAAAGAAA	92379	92419	1	-	0	target=ACAAGAATATTTTTTTTATTGGGACGATAATTCTGAAAAGAAAGA
chrc	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:12568:1114#0/1:ATAAATTAGAATTAAAATTTTAATAGTTGTAATTAATATTG	112309	112349	1	+	0	target=CTATAAATCAGAACCAAAATCCCAACAGTTGTAATTAATATTGAC
.	NM/NM	SOLEXA2_0531_FC62W31AAXX:4:1:12632:1114#0/1:ATGTGTATGATTGAGTATAAGAATTTAAATCGGAATTATAT	0	0	0	.	.	.
chr2	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:12896:1114#0/1:GTAATGGATGTTTGAGGTAAAAATGTTTTATTTTTTTTTGT	11592020	11592060	1	+	0	target=CTGTAATGGATGCCTGAGGCAAAAATGCTCCATTTTCTTTTGCGC
.	NM/NM	SOLEXA2_0531_FC62W31AAXX:4:1:12999:1112#0/1:AATTAATTATTTTTTTGTTTTTTNAAGTTTTTATGGTGTAG	0	0	0	.	.	.
chr1	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:13022:1111#0/1:TGTTGTTGATGTTTATTTATAATNTAAATAAAATTAAGTTT	12903623	12903663	1	-	1	target=GTTGTTGTTGATGTTTATCTATAATACAAACAAAATTAAGTTCAA
chr4	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:13303:1114#0/1:TAATTTTTTGTGAGATTGATTTTTTTTGATATTGATGAATT	7807569	7807609	1	+	0	target=AACAATCTTCCGTGAGATCGATTTCTTCTGACATCGATGAATCGA
chr4	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:13332:1112#0/1:GTAGTATGTTTTTAGGAAAAGTGNTTATGTATTTTGGTTTA	5181385	5181425	1	-	1	target=AAGTAGTATGCCTTCAGGAAAAGCGCTCATGTATTTCGGTTTACT
chr1	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:13432:1114#0/1:GATGGATAAGGATGATTGTTTGGATGTTAAAGAGATTGTTT	26335637	26335677	1	+	0	target=GCGACGGATAAGGATGATCGTCCGGATGCTAAAGAGATTGTTCAG
.	NM/NM	SOLEXA2_0531_FC62W31AAXX:4:1:13459:1117#0/1:TTTTTGGTTTTGTGTTTGTTAATAAGGATATATTATTTAGG	0	0	0	.	.	.
chr3	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:13594:1115#0/1:GAGATTTTTTTTTAAAAGTATTTAGTTGGTGGGATTTTTTT	8640238	8640278	1	-	0	target=TGGAGATTCTCCTTCAAAAGCATTCAGTTGGTGGGACTCTTTCCC
chr5	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:13641:1117#0/1:TATTTGTTTAATATTAGTATTTGTTAAGAAATATATTTTAG	11440858	11440898	1	-	0	target=ATTATCTGTTTAATATTAGTATTTGTTAAGAAATATATTTTAGTT
chr5	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:13713:1110#0/1:AAAATTTTATGATTCTACAATTTNTAAACAGTGAGGTCAAA	14110325	14110365	1	-	1	target=TTAAAACTTTATGACTCTACAACTTTTAAACAGTGAGGCCAAACT
chrc	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:13807:1112#0/1:TTTTATTTAATTTTAAATTAATGAGTTATATATTAAGTATA	82613	82653	1	+	0	target=ATCTTTATCTAATCCTAAACCAACGAGTCACACACTAAGCATAGC
chrm	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:13863:1110#0/1:AGAATTTGATTTTGTATTTTTGTNTAATTGATATATGATAT	321898	321938	1	+	1	target=AAAGAATTCGACCTTGCACTCCCGCTCAATCGATATACGACATCG
chr5	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:13896:1114#0/1:TAATTTTATTGTTTGTTAGTTTTTTTTGTGAATTTTTTTGG	2145097	2145137	1	+	0	target=TTCAACTTCATTGTCTGTTAGTCTCTCTTGTGAATTTCTTTGGAG
chr1	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:14156:1115#0/1:TTTTAAGATTTGAGTTCTGTTATTGATGATTATTTCGATAG	8367571	8367611	1	-	0	target=ATCTTCAAGATTTGAGTTCTGTTATTGATGACTATCTCGATAGCC
chr1	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:14197:1113#0/1:TTTTTGAAGAAATTAGTGATTAAAGAATATTTTTTGAATTT	713775	713815	1	-	0	target=GATTCTTGAAGAAATCAGTGACCAAAGAACATTTCTTGAATTTGA
chrc	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:14261:1109#0/1:GTATTTAAGTAGTAAGTTTATTTNAAGATGAGTGTTTTTTT	107460	107500	1	+	1	target=AAGCATCTAAGTAGTAAGCCCACCCCAAGATGAGTGCTCTCCTAT
chrm	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:14445:1110#0/1:TAATTTTTGGAATTGTTAATTTGNATTTAGTTAATAGTTGA	196394	196434	1	+	2	target=GCTAACTCCCGGAATCGCCAACCCGAATTCAGCTAACAGCTGCTA
chrm	R/NM*	SOLEXA2_0531_FC62W31AAXX:4:1:14606:1116#0/1:TTTATATTTTTAATTTGTTGGTTTGTTAGGAAGGGAAGTTA	118215	118255	1	-	0	target=TCTTTACATCTCTAACTTGTTGGCTCGTTAGGAAGGGAAGTTAAA
chrc	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:14699:1111#0/1:GTTATGAATATAAAGTTTTAAAGNATGGAATTTTAGAAAGA	38260	38300	1	+	1	target=TCGTTATGAACATAAAGTCCCAAAGTATGGAACCCTAGAAAGAGG
chr5	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:14896:1111#0/1:GATTTGTGTTTTGAATGGTTTTGNTAAAACGAATAATATTG	19981785	19981825	1	-	1	target=GAGACCTGCGTCTTGAATGGTTTTGATAAAACGAACAACATCGAG
chr4	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:14919:1110#0/1:AGGGTTTATTTGTAATTAATTTTNTAAATTTGTTTGGTATG	6865041	6865081	1	+	1	target=TAAGGGTTTACTTGTAACCAATCCTTCAAACCCGCTTGGTACGAC
chrc	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:15101:1115#0/1:GAAAAAGTTGTATTTTTGATGGAAGAATAGGAGATTTTTTT	130923	130963	1	-	0	target=AGGAAAAAGCCGCATTTTTGATGGAAGAACAGGAGATCCTTTTGA
chr4	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:15374:1111#0/1:ATTCGAAAATGATGAAATTCTAGNATGAAAATGCAACAAGA	4731819	4731859	1	+	1	target=AAATTCGAAAATGATGAAACCCTAGCATGAAAATGCAACAAGAAA
chrc	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:15392:1114#0/1:TGTTTTTATTATGGATGTTTTATTTTGTTTTTATTTATGAA	75664	75704	1	-	0	target=AGCGTCTTTATTATGGACGCTTTATTCTGTCTCCACTTATGAAAG
chr4	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:15848:1115#0/1:TGTAGTTTTGATTAGATTTTGGTAATGTAATTTTGAGAAAT	13102533	13102573	1	+	0	target=TTTGTAGCTCTGACCAGACTTTGGCAACGTAATTTTGAGAAATGC
.	NM/NM	SOLEXA2_0531_FC62W31AAXX:4:1:15967:1116#0/1:GAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATG	0	0	0	.	.	.
chr1	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:16082:1113#0/1:TTTTATGTGATTTTAATAAAGTTATAGTTTTAAAAGTGGAA	8517878	8517918	1	+	0	target=GTTTCTATGTGACTTTAATAAAGTTATAGTTTCAAAAGTGGAAGC
chrc	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:16250:1112#0/1:TATAAGTGGTTGAATTTTTATTANGTAAGGGTTTATTGGAT	58417	58457	1	+	1	target=TTCACAAGCGGCTGAATCTTTATTACGTAAGGGCTTATTGGATGC
chr5	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:16407:1112#0/1:TGTGTTTGTTTAAATTTTTATATTATTTTGGTATTTTTGTA	13912450	13912490	1	+	0	target=TTTGTGTTTGTTCAAACTCTTACATTACTTTGGTATCTTTGTATA
.	NM/NM	SOLEXA2_0531_FC62W31AAXX:4:1:16865:1112#0/1:GAAGAGCGGTTCAGCAGGAATGCNGAGACCGATCTCGTATG	0	0	0	.	.	.
chr2	U/NM	SOLEXA2_0531_FC62W31AAXX:4:1:16907:1112#0/1:ATTTTTTGTAGTGTTAAAATAAANATAGGTGAATTTATTGA	13283614	13283654	1	+	1	target=CTATTCTTCGCAGTGTTAAAACAAAAATAGGTGAATTTACTGAAA
.	NM/NM	SOLEXA2_0531_FC62W31AAXX:4:1:17033:1114#0/1:GTGTTTTTTTGTTAGAAGACATAAAGTTAAAGATTTATATG	0	0	0	.	.	.
