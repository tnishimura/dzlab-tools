#!/usr/bin/env perl
use Getopt::Long;
use strict;
use warnings;
#use diagnostics;disable diagnostics;
use Carp;
use Data::Dumper;
use List::Util qw(min);
use feature 'say';
use autodie;

use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if (!$opt_reference || !$opt_eland);


# holds name of chromosomes as keys and length of chromosomes in bp as values
my %reference = %{ index_fasta ($opt_reference) };

# redirects STDOUT to file if specified by user
if ( $opt_output ne q{-} ) {
    open STDOUT, '>', "$opt_output";
}

# opens sequence files
open my $LEFT,  '<', $opt_eland
or croak "Can't open file: $opt_eland";


while (my $leftend = <$LEFT>) {
    # reads single sequence from each file

    $leftend =~ s/[\r\n]//g;

    # parses each line into hash
    my %left  = %{ parseEland3Line($leftend, $opt_max_hits) };

    # no matches
    # $VAR1 = {
    #     'sequence' => 'GTTTTGGAGTCGAATATGATTTGATGTTATGTGTATGATTGAGTA',
    #     'line' => 'HWI-EAS105:7:1:7:715#0/1'
    #     'matches' => 0,
    # };

    # 1 match
    # $VAR1 = {
    #     'sequence' => 'GTTTTAAGTTGTTGGGTTGTTNAGTGTTTAAATATTTGATTAAAT',
    #     'line' => 'HWI-EAS105:7:1:7:1958#0/1'
    #     'matches' => 1,
    #     'coord0' => '29272572',
    #     'chr0' => 'RC_Chr1',
    #     'mm0' => '1',
    #     'rank' => [ '0', '1', '0', '0' ],
    # };

    # > 1 matches
    # $VAR1 = {
    #     'sequence' => 'TTTAATATTGATGAAATTAAANTAGGAATTTTTTGGGAAAAAGGT',
    #     'line' => 'HWI-EAS105:7:1:7:565#0/1'
    #     'matches' => 2,
    #     'coord0' => '3329121',
    #     'mm0' => '1',
    #     'chr0' => 'Chr2',
    #     'mm1' => '1',
    #     'chr1' => 'RC_ChrM',
    #     'coord1' => '29638',
    #     'rank' => [ '0', '2', '0', '0' ],
    # };

    # gets number of matches and sequence from hash
    my $lmatch    = $left{'matches'};
    my $lsequence = $left{'sequence'};
    my $lreadsize = length $lsequence;

    # initializes each gff field to default values
    my (
        $l_seqname, $l_source, $l_feature, $l_start, $l_end,
        $l_score,   $l_strand, $l_frame,   $l_attribute
    ) = ( q{.}, 'pair_align', q{.}, 0, 0, 0, q{.}, q{.}, q{.} );

    if ( $lmatch == 0) {
        $l_feature = $left{'line'} . q{:} . $lsequence;
        $l_source  = 'NM/NM';
    }

    if ( $lmatch == 1) {

        $l_source  = 'U/NM';

        # gets chromosome name into $tmp
        $left{'chr0'} =~ m/(.*)/i;
        my $tmp = $1;
        $tmp =~ tr/A-Z/a-z/;
        $tmp =~ s/rc_//i;

        if ( $left{'chr0'} =~ m/^RC_/i ) {  # if left sequence maps to reverse strand
            $l_seqname = $tmp;
            $l_feature = $left{'line'} . q{:} . $lsequence;
            $l_start   = $left{'coord0'};
            $l_end     = $left{'coord0'} + $lreadsize - 1;
            $l_score   = 1;
            $l_strand  = q{-};
            $l_frame   = $left{mm0};
            $l_attribute =
            'target='
            . substr( $reference{"$tmp-rc"}, $left{'coord0'} - 1,
                      $lreadsize + 4 );
        }
        else { # if left sequence maps to forward strand

            $l_seqname   = $tmp;
            $l_feature   = $left{'line'} . q{:} . $lsequence;
            $l_start     = $left{'coord0'};
            $l_end       = $left{'coord0'} + $lreadsize - 1;
            $l_score     = 1;
            $l_strand    = q{+};
            $l_frame     = $left{mm0};
            $l_attribute = 'target='
            . substr(
                $reference{"$tmp-seq"},
                $left{'coord0'} - 1,
                $lreadsize + 4
            );
        }
    }

    if ( $lmatch > 1) {
        $l_source = "R/NM";

        $l_feature   = $left{'line'} . q{:} . $lsequence;
        $l_score     = $lmatch;

        if ($opt_random_assign) {

            $l_source .= q{*};

            my $rnd = int rand $lmatch;
            
            # gets chromosome name into $tmp
            $left{"chr$rnd"} =~ m/(.*)/i;
            my $tmp = $1;
            $tmp =~ tr/A-Z/a-z/;
            $tmp =~ s/rc_//i;

            if ( $left{"chr$rnd"} =~ m/^RC_/i ) { # if left sequence maps to reverse strand
                $l_seqname   = $tmp;
                $l_feature   = $left{"line"} . q{:} . $lsequence;
                $l_start     = $left{"coord$rnd"};
                $l_end       = $left{"coord$rnd"} + $lreadsize - 1;
                $l_score     = 1;
                $l_strand    = q{-};
                $l_frame     = $left{"mm$rnd"};
                $l_attribute = 'target='
                . substr(
                    $reference{"$tmp-rc"},
                    $left{"coord$rnd"} - 1,
                    $lreadsize + 4
                );
            }
            else {        # if left sequence maps to forward strand
                $l_seqname   = $tmp;
                $l_feature   = $left{"line"} . q{:} . $lsequence;
                $l_start     = $left{"coord$rnd"};
                $l_end       = $left{"coord$rnd"} + $lreadsize - 1;
                $l_score     = 1;
                $l_strand    = q{+};
                $l_frame     = $left{"mm$rnd"};
                $l_attribute = 'target='
                . substr(
                    $reference{"$tmp-seq"},
                    $left{"coord$rnd"} - 1,
                    $lreadsize + 4
                );
            }
        }
        else {

            $l_attribute = 'targets=';

            # gets all possible matches and appends them to $r_attribute
          NM_MM_LOOP_LEFT:
            for my $i ( 0 .. $lmatch - 1 ) {
                $left{"chr$i"} =~ m/(.*)/i;
                my $tmp = $1;
                $tmp =~ tr/A-Z/a-z/;
                $tmp =~ s/rc_//i;

                if ( $left{"chr$i"} =~ m/^RC_/i ) { # if left sequence maps to reverse strand

                    # appends chromosome name:coordinate:sequence to $l_attribute
                    $l_attribute .= join(
                        q{:},
                        $left{"chr$i"},
                        $left{"coord$i"},
                        substr(
                            $reference{"$tmp-rc"},
                            $left{"coord$i"} - 1,
                            $lreadsize + 4
                        )
                    );

                } else {     # if left sequence maps to reverse strand

                    # appends chromosome name:coordinate:sequence to $l_attribute
                    $l_attribute .= join(
                        q{:},
                        $left{"chr$i"},
                        $left{"coord$i"},
                        substr(
                            $reference{"$tmp-seq"},
                            $left{"coord$i"} - 1,
                            $lreadsize + 4
                        )
                    );
                }
                $l_attribute .= q{,}; # appends sub-field separator comma
            }
            chop $l_attribute; # deletes very last comma and appends new line to $r_attribute
        }
    }

    # print both ends
    print join( "\t",
                $l_seqname, $l_source, $l_feature, $l_start, $l_end,
                $l_score,   $l_strand, $l_frame,  $l_attribute ),
                "\n"
                unless ($opt_skip_nm and $l_score = 1);
}

# end for loop through all sequences
print STDERR "done!\n";
close STDOUT;
close $LEFT;

# main program finished

# HWI-EAS105:7:1:7:1958#0/1	GTTTTAAGTTGTTGGGTTGTTNAGTGTTTAAATATTTGATTAAAT	0:1:0:0	RC_Chr1:29272572F1
# HWI-EAS105:7:1:7:565#0/1	TTTAATATTGATGAAATTAAANTAGGAATTTTTTGGGAAAAAGGT	0:2:0:0	Chr2:3329121F1,RC_ChrM:29638F1
# HWI-EAS105:7:1:7:706#0/1	ACGGTTATGCGTATTTTTTAANTATTTGTGTCAGAGAGTTATTGA	0:1:0:0	RC_Chr2:15907631F1
# HWI-EAS105:7:1:7:221#0/1	GAATTTGAAAACGATAATAATNTTAAAATTTTGGAGGAAGGATGG	0:2:0:0	Chr5:15705316F1,Chr5:15691740F1
# HWI-EAS105:7:1:7:1851#0/1	TAAAATTTTAGATTTGAGTTGNATTTTTTTGATGTTTGGAATATT	0:1:0:0	Chr4:12792490F1
# HWI-EAS105:7:1:7:1125#0/1	GTAAAATAGGTTTTTTTGAGANAAATTTATGAGAAATATAATGGT	0:1:0:0	Chr1:2643593F1
# HWI-EAS105:7:1:7:876#0/1	TTGAAGATTGTTATAATTTAANATTGAATTTAATGTTAAGAAATT	0:1:0:0	Chr3:23270809F1
# HWI-EAS105:7:1:7:1859#0/1	GTTAGTTGAAGAGTAAATCGAAGAAAATATGAGGAAGGGTTATGA	1:0:0:0	RC_Chr2:17917881F0
# HWI-EAS105:7:1:7:1385#0/1	TTTATGATTAAGGAAGAGTTTGTTAAGGGAGTTAAGTTGTTGTTT	1:0:0:0	RC_Chr2:3014191F0
# HWI-EAS105:7:1:7:1402#0/1	GATTAATTAGGATGTGTTTATGGATGAAATAATTTTAGATTTTGA	1:0:0:0	RC_Chr3:584697F0

# Parses single eland3 record (output from seqmap)
# Returns reference to hash with following keys:
# 'line', 'sequence', 'matches', 'chr[0-#matches-1]', 'coord[0-#matches-1]'
sub parseEland3Line {
    my ($tmp_eland_line, $max_hits) = @_;
    $tmp_eland_line =~ s/[\n\r]//g;
    my @eland_line = (split /\t/, $tmp_eland_line);

    my %hash = ();
    $hash{'line'}     = $eland_line[0];
    $hash{'sequence'} = $eland_line[1];

    if ($eland_line[2] =~ m/NM/) {
        $hash{'matches'} = 0;
    }

    elsif ($eland_line[2] =~ m/^[0-9]+$/) {
        $hash{'matches'} = $eland_line[2];
        $hash{rank}      = $eland_line[2];
    }

    elsif ($eland_line[2] =~ m/:/) {
        my @all_matches   = split /,/, $eland_line[3];
        $hash{'matches'}  = scalar @all_matches;
        @{$hash{'rank'}}  = split /:/, $eland_line[2];
    }

    if ($hash{'matches'} > 1) {

        my @all_reads
        = sort { (substr $a, -1) <=> (substr $b, -1) }
        (split /,/, $eland_line[3]);

        if ($max_hits) {

          RANK:
            for my $rank (0 .. @{$hash{rank}} - 1) {
                
                if ($hash{rank}->[$rank] == 0) {
                    $hash{matches} -= $hash{rank}->[$rank];
                    next RANK;
                }
                elsif ($hash{rank}->[$rank] > $max_hits) {
                    $hash{matches} = 0;
                    return \%hash;
                }

              RANKHIT:
                for my $i (0..@all_reads - 1) {

                    my ($tmp_chr, $tmp_coord) = split /:/, $all_reads[$i];
                    $tmp_coord =~ s/[A-Z]([0-9])$//i;
                    my $tmp_mm = $1;

                    if ($tmp_mm == $rank) {
                        ($hash{"chr$i"}, $hash{"coord$i"}, $hash{"mm$i"})
                        = ($tmp_chr, $tmp_coord, $tmp_mm);
                    }
                }
            }
        }
        else {
          HIT:
            for my $i (0 .. @all_reads - 1) {
                ($hash{"chr$i"}, $hash{"coord$i"}) = split /:/, $all_reads[$i];
                $hash{"coord$i"} =~ s/[A-Z]([0-9])$//i;
                $hash{"mm$i"} = $1;
            }
        }
    }
    elsif ($hash{'matches'} == 1) {
        ($hash{'chr0'}, $hash{'coord0'}) = split /:/, $eland_line[3];
        $hash{'coord0'} =~ s/[A-Z]([0-9])$//i;
        $hash{'mm0'} = $1;
    }
    return \%hash;
}



# Converts input sequence to reverse complement
# Returns scalar $string with processed sequence
sub reverseComp {
    my $shortseq = shift;
    $shortseq =~ tr/ACGTacgt/TGCAtgca/;
    $shortseq =~ s/\n//;
    return reverse $shortseq;
}


sub index_fasta {
    my ($referencefile) = @_;

    # holds name of chromosomes as keys and length of chromosomes in bp as values
    my %reference = ();

    # reads in the reference genome file into @fastaseq
    open my $REF, '<', "$referencefile" or croak "Can't open file: $referencefile";
    my @fastaseq = <$REF>;
    close $REF;

    # find and store indices for each chromosome change and corresponding descriptions
    my ( @idx, @dsc ) = ();
    for my $i ( 0 .. @fastaseq - 1 ) { ### Indexing $referencefile...  % done
        if ( $fastaseq[$i] =~ m/^>/ ) {
            $fastaseq[$i] =~ s/>//g;
            $fastaseq[$i] = ( split /\s/, "$fastaseq[$i]" )[0];
            $fastaseq[$i] =~ tr/A-Z/a-z/;
            push @idx, $i;
            push @dsc, $fastaseq[$i];
        }
    }

    # gets and saves each chromosome's sequence and reverse complemented sequence
    for my $j ( 0 .. @idx - 1 ) { ### Loading $referencefile into memory...  % done
        my $line;
        if ( $j == scalar @idx - 1 ) {
            $line = join( q{}, 'NN', @fastaseq[ $idx[$j] + 1 .. @fastaseq - 1] , 'NN');
        }
        else {
            $line = join( q{}, 'NN', @fastaseq[ $idx[$j] + 1 .. $idx[$j + 1] - 1] , 'NN');
        }
        $line =~ s/[\n\r]//g;
        $reference{ $dsc[$j] }     = (length $line) - 4;
        $reference{"$dsc[$j]-seq"} = $line;
        $reference{"$dsc[$j]-rc"}  = reverseComp ($line);
    }
    return \%reference;
}


=head1 NAME

correlateSingleEnds.pl - ...

=head1 SYNOPSIS

Usage examples:

 correlateSingleEnds.pl [options]...

=head1 REQUIRED ARGUMENTS

=over

=back

=head1 OPTIONS

=over

=item  -e <eland> | --eland <eland>

Output of parse_bowtie.pl

=for Euclid
    eland.type:        readable

=item  -f <fasta> | --reference <fasta>

Genome reference.

=for Euclid
    fasta.type:        readable

=item  -rnd <rnd> | --random-assign <rnd>

If a read matches to multiple places, pick one random. 

=for Euclid
    rnd.default:     1

=item --help | -h

=item  -o <file> | --output <file>

Output file

=for Euclid
    file.default:     '-'

=item  -k  | --skip-nm 

=item  -m <hits> | --max-hits <hits>

=for Euclid
    hits.default:     0

=back

=cut

__DATA__

HWI-EAS105:7:1:7:1958#0/1	GTTTTAAGTTGTTGGGTTGTTNAGTGTTTAAATATTTGATTAAAT	0:1:0:0	RC_Chr1:29272572F1
HWI-EAS105:7:1:7:565#0/1	TTTAATATTGATGAAATTAAANTAGGAATTTTTTGGGAAAAAGGT	0:2:0:0	Chr2:3329121F1,RC_ChrM:29638F1
HWI-EAS105:7:1:7:706#0/1	ACGGTTATGCGTATTTTTTAANTATTTGTGTCAGAGAGTTATTGA	0:1:0:0	RC_Chr2:15907631F1
HWI-EAS105:7:1:7:221#0/1	GAATTTGAAAACGATAATAATNTTAAAATTTTGGAGGAAGGATGG	0:2:0:0	Chr5:15705316F1,Chr5:15691740F1
HWI-EAS105:7:1:7:1851#0/1	TAAAATTTTAGATTTGAGTTGNATTTTTTTGATGTTTGGAATATT	0:1:0:0	Chr4:12792490F1
HWI-EAS105:7:1:7:1125#0/1	GTAAAATAGGTTTTTTTGAGANAAATTTATGAGAAATATAATGGT	0:1:0:0	Chr1:2643593F1
HWI-EAS105:7:1:7:876#0/1	TTGAAGATTGTTATAATTTAANATTGAATTTAATGTTAAGAAATT	0:1:0:0	Chr3:23270809F1
HWI-EAS105:7:1:7:1859#0/1	GTTAGTTGAAGAGTAAATCGAAGAAAATATGAGGAAGGGTTATGA	1:0:0:0	RC_Chr2:17917881F0
HWI-EAS105:7:1:7:1385#0/1	TTTATGATTAAGGAAGAGTTTGTTAAGGGAGTTAAGTTGTTGTTT	1:0:0:0	RC_Chr2:3014191F0
HWI-EAS105:7:1:7:1402#0/1	GATTAATTAGGATGTGTTTATGGATGAAATAATTTTAGATTTTGA	1:0:0:0	RC_Chr3:584697F0

chr1	U/U	HWI-EAS105:7:1:7:1958#0/1:GTTTTAAGTTGTTGGGTTGTTNAGTGTTTAAATATTTGATTAAAT	29272572	29272616	1	-	1	target=AGGCTTTAAGTTGTTGGGCTGTTTAGTGTTCAAACATTTGATTAAACCG
chr1	U/U	HWI-EAS105:7:1:7:1958#0/1:TGATTAAAAAGTTATTATATATTAATAGTAG	29272617	29272647	1	-	2	target=ACCGATTAAAAAGCTACTATACATTAACAGCCAAG
chr2	R/R	HWI-EAS105:7:1:7:565#0/1:TTTAATATTGATGAAATTAAANTAGGAATTTTTTGGGAAAAAGGT	3329121	3329165	1	+	1	target=TATTCAATATCGACGAAATTAAATTAGGAATCTTTTGGGAAAAAGGCCC
chr2	R/R	HWI-EAS105:7:1:7:565#0/1:TTTTTTTAGAGGTTGTTAATTTATGAATTAG	3329166	3329196	1	+	0	target=GCCCTTTCTAGAGGTTGTCAATTTACGAATTAGTG
chr2	U/U	HWI-EAS105:7:1:7:706#0/1:ACGGTTATGCGTATTTTTTAANTATTTGTGTCAGAGAGTTATTGA	15907631	15907675	1	-	1	target=CTACGGCTATGCGCATTCTTCAATTATTTGTGTCAGAGAGTTATTGACA
chr2	U/U	HWI-EAS105:7:1:7:706#0/1:TATTGAGGTTATTTTAGTATTGGAAATTGAA	15907676	15907706	1	-	0	target=GACATTGAGGTTATTTCAGTATTGGAAACTGAAAT
chr5	R/NM*	HWI-EAS105:7:1:7:221#0/1:GAATTTGAAAACGATAATAATNTTAAAATTTTGGAGGAAGGATGG	15705316	15705360	1	+	1	target=TAGAATCTGAAAACGATAACAACTTTAAAACTTTGGAGGAAGGATGGCC
.	R/NM	HWI-EAS105:7:1:7:221#0/1:CCGAGACGATAATGTTAAGATCGGAAGAGCG	0	0	0	.	.	.
chr4	U/U	HWI-EAS105:7:1:7:1851#0/1:TAAAATTTTAGATTTGAGTTGNATTTTTTTGATGTTTGGAATATT	12792490	12792534	1	+	1	target=TCCAAAACCTCAGACTTGAGCTGGACTTTCTCGATGTCTGGAACATCTG
chr4	U/U	HWI-EAS105:7:1:7:1851#0/1:TGGATGAGTGAGGAGGAAGAGTTTGTTGGTT	12792535	12792565	1	+	0	target=TCTGGATGAGTGAGGAGGAAGAGCTTGTTGGCTAG
chr1	U/U	HWI-EAS105:7:1:7:1125#0/1:GTAAAATAGGTTTTTTTGAGANAAATTTATGAGAAATATAATGGT	2643593	2643637	1	+	1	target=CCGTAAAACAGGTTCTCCTGAGAAAAATTTACGAGAAACATAATGGTCA
chr1	U/U	HWI-EAS105:7:1:7:1125#0/1:TAAGTTTTTTTGTGAGATTTTTTTATAAGTT	2643638	2643668	1	+	0	target=GTCAAGCCCTTCCGCGAGATTTCTCCACAAGTCCA
chr3	U/U	HWI-EAS105:7:1:7:876#0/1:TTGAAGATTGTTATAATTTAANATTGAATTTAATGTTAAGAAATT	23270809	23270853	1	+	1	target=TTCTGAAGACTGTTACAACCTAAAATTGAATTCAATGTTAAGAAACCCC
chr3	U/U	HWI-EAS105:7:1:7:876#0/1:TTGAAAAAAGTTGTAAATTTTGAATATTAAG	23270854	23270884	1	+	1	target=CCCCGAAAAAAGCTGCAAATTTTGAACACTAAAAG
chr2	U/U	HWI-EAS105:7:1:7:1859#0/1:GTTAGTTGAAGAGTAAATCGAAGAAAATATGAGGAAGGGTTATGA	17917881	17917925	1	-	0	target=AAGTCAGTTGAAGAGTAAACCGAAGAAAACATGAGGAAGGGTCACGAGA
chr2	U/U	HWI-EAS105:7:1:7:1859#0/1:GAGATTTGAGAGGGAGATATGGAGAAAGATT	17917926	17917956	1	-	0	target=GAGAGACTTGAGAGGGAGACATGGAGAAAGATCGC
chr2	U/NM	HWI-EAS105:7:1:7:1385#0/1:TTTATGATTAAGGAAGAGTTTGTTAAGGGAGTTAAGTTGTTGTTT	3014191	3014235	1	-	0	target=TGTTCATGATCAAGGAAGAGCTTGTCAAGGGAGCTAAGCTGTTGCTCTC
.	U/NM	HWI-EAS105:7:1:7:1385#0/1:TTAAATGTTTAGTTTTTTTGATAAGTTTTTT	0	0	0	.	.	.
chr3	U/NM	HWI-EAS105:7:1:7:1402#0/1:GATTAATTAGGATGTGTTTATGGATGAAATAATTTTAGATTTTGA	584697	584741	1	-	0	target=GGGACCAACCAGGATGTGTCTATGGATGAAACAACTTCAGATTCTGATC
.	U/NM	HWI-EAS105:7:1:7:1402#0/1:TTTTAATAGATCGGAAGAGCGGTTCAGCAGG	0	0	0	.	.	.
