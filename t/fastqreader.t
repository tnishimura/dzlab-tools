#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More qw(no_plan);
use Test::Exception;
use DZUtil qw/c2t g2a reverse_complement/;
use FastqReader;
use FastqReader::Convert;
use FastqReader::CountReads;
use FastqReader::GetReads;
use FastqReader::IsBS;
my $count = 5100;
my $fq   = "t/data/bs-sequel-test.fastq";

my $fa     = "t_intermediate/fasta.fasta";
my $fa_c2t = "t_intermediate/c2t.fasta";
my $fa_g2a = "t_intermediate/g2a.fasta";
my $fa_rc     = "t_intermediate/rc.fasta";
my $fa_c2t_rc = "t_intermediate/c2t.rc.fasta";
my $fa_g2a_rc = "t_intermediate/g2a.rc.fasta";

my $fq_c2t = "t_intermediate/c2t.fastq";
my $fq_g2a = "t_intermediate/g2a.fastq";
my $fq_c2t_rc = "t_intermediate/c2t.rc.fastq";
my $fq_g2a_rc = "t_intermediate/g2a.rc.fastq";

fastq_convert(in => $fq, out => $fa);
fastq_convert(methyl => 'c2t',                         in => $fq, out => $fa_c2t);
fastq_convert(methyl => 'g2a',                         in => $fq, out => $fa_g2a);
fastq_convert(methyl => 'c2t', to_fasta => 0,          in => $fq, out => $fq_c2t);
fastq_convert(methyl => 'g2a', to_fasta => 0,          in => $fq, out => $fq_g2a);
fastq_convert(methyl => 'c2t',                rc => 1, in => $fq, out => $fa_c2t_rc);
fastq_convert(methyl => 'g2a',                rc => 1, in => $fq, out => $fa_g2a_rc);
fastq_convert(methyl => 'c2t', to_fasta => 0, rc => 1, in => $fq, out => $fq_c2t_rc);
fastq_convert(methyl => 'g2a', to_fasta => 0, rc => 1, in => $fq, out => $fq_g2a_rc);

my $fqr = FastqReader->new(file => \*DATA, fasta => 0);

is(count_reads($fq),           $count, "read count of fq"); 
is(count_reads($fq_c2t),       $count, "read count of fq c2t"); 
is(count_reads($fq_g2a),       $count, "read count of fq g2a"); 
is(count_reads_fasta($fa),     $count, "read count of fasta"); 
is(count_reads_fasta($fa_c2t), $count, "read count of fasta c2t"); 
is(count_reads_fasta($fa_g2a), $count, "read count of fasta g2a"); 

check_conversion(0, 1, 0, $fq, $fa_c2t);
check_conversion(0, 0, 1, $fq, $fa_g2a);
check_conversion(0, 1, 0, $fq, $fq_c2t);
check_conversion(0, 0, 1, $fq, $fq_g2a);
check_conversion(1, 1, 0, $fq, $fa_c2t_rc);
check_conversion(1, 0, 1, $fq, $fa_g2a_rc);
check_conversion(1, 1, 0, $fq, $fq_c2t_rc);
check_conversion(1, 0, 1, $fq, $fq_g2a_rc);

no_c($fq_c2t);
no_c($fq_c2t_rc);

sub no_c{
    my $fastq = shift;
    my ($is_bs, %ratios) = is_bs($fastq);
    ok($is_bs, "$fastq is_bs");
    is(0, $ratios{c}, "$fastq c_ratio == 0");
}

sub check_conversion{
    my ($rc, $c2t, $g2a, $original_file, $converted_file) = @_;
    my $fqr_original  = FastqReader->new(file => $original_file,  fasta => scalar($original_file =~ /fasta$/));
    my $fqr_converted = FastqReader->new(file => $converted_file, fasta => scalar($converted_file =~ /fasta$/));

    my $line = 1;
    my @errors;
    
    while (!$fqr_original->eof() && !$fqr_converted->eof()){
        my $o = $fqr_original->next();
        my $c = $fqr_converted->next();

        my $seq_o = $o->[1];
        my $seq_c = $c->[1];

        $seq_o = reverse_complement($seq_o) if $rc;
        $seq_o = c2t($seq_o) if $c2t;
        $seq_o = g2a($seq_o) if $g2a;

        if ($seq_o ne $seq_c){
            push @errors, $line;
        }

        $line += 4;
    }
    ok(scalar(@errors) == 0, "$original_file => $converted_file, rc = $rc, c2t = $c2t, g2a = $g2a");
}

#######################################################################

{
    my $reads = get_reads(\*DATA, "HS2:90:B09PCABXX:1:1202:13269:194742",);
    is(scalar(@$reads), 1, "get_reads len = 1");
    is($reads->[0][1], "ATTTTTTTATTTTTGTGATCCGGTTGCGGTTTAAGTTGTTATATTTAATGATATACAGGATATTAAGTTATATTTGATTTTAAAAATTTAATTAATTTTT", "read seq");
}

seek \*DATA, 0,0;

{
    my $reads = get_reads(\*DATA, '@HS2:90:B09PCABXX:1:1202:13292:194744#/1 1:@:3',);
    is(scalar(@$reads), 1, "get_reads len = 1");
    is($reads->[0][1], "GTTTAGATGTAGTGGTTGTGAAAGTTAAAATGTAGAAGGTTGATAATTTATTTGAATTAATTTGATTTTGTTGGAAGTTTGATATTTTTGGAAAAAAGTG", "read seq");
}

__DATA__
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
