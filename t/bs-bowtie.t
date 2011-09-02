#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Launch;
use File::Basename;
use File::Spec::Functions;
use Test::More tests => 14;
use TestUtils;
use FastaReader;
use DZUtil qw/numdiff c2t/;

my $test_dir = 't_intermediate';
my $ref = setup_reference($test_dir);

my $num_reads = 10000;
my $num_junk = 1000;
my $total = $num_reads + $num_junk;
my $bowtieout = catfile($test_dir, "bs-bowtie-test");
my $bowtieout_left = catfile($test_dir, "bs-bowtie-test_1-60");
my $bowtieout_right = catfile($test_dir, "bs-bowtie-test_61-100");
my $testfastq = catfile($test_dir, "bs-bowtie-sheared.fastq");
my $refreader = FastaReader->new(file => $ref, slurp => 1);

ok(launch("perl ./genome_shear.pl -j $num_junk -n $num_reads -l 100 -o $testfastq $ref"),
    'genome shear ran');
ok(launch("perl ./bs-bowtie -f $ref -r $testfastq -o $bowtieout"), 'whole read bs-bowtie ran');
ok(launch("perl ./bs-bowtie -s 1 60 -f $ref -r $testfastq -o $bowtieout_left"), 'left read bs-bowtie ran');
ok(launch("perl ./bs-bowtie -s 61 100 -f $ref -r $testfastq -o $bowtieout_right"), 'right read bs-bowtie ran');

for my $bout ($bowtieout, $bowtieout_left, $bowtieout_right){
    like(`wc -l $bout`, qr/$total/, "correct number of reads output from bs-bowtie");
    ok(junk_count($bout) >= $num_junk, "bowtie output has at least as many junks as inputed");
}
ok(reads_match(100, $bowtieout, $refreader), 'whole reads match alignment location');
ok(reads_match(60, $bowtieout_left, $refreader), 'left reads match alignment location');
ok(reads_match(40, $bowtieout_right, $refreader), 'right reads match alignment location');

# chr3:563482:563581:-:563518:563519:563573#/1	GATAAGTTCTAAAAAAGTAAAATTGATAAATGGTTTTAATTTATTTTTTTTTTGTTTTTT	1:0:0:0	RC_chr3:1016341F0
# chr5:187766:187865:+:187768:187771:187796:187835#/1	ATCGGCAGTGTTTAGTATATATAGTTAATGCTGTGTTTGAAAAGATATGGTTTATGTGTA	1:0:0:0	chr5:187766F0
# chr4:111873:111972:-:111874#/1	TTAAGAAAGTTATAGAATTGTATTATGATATTGTGAAGTTGTAGAAAAATAAATAATAGT	1:0:0:0	RC_chr4:1467950F0
# chr4:457333:457432:-:457385#/1	ATATATTTAAATTGAATTTTGATGAATGAAAAGATGGGAAGAATTTGCTTTAGTAAATAG	1:0:0:0	RC_chr4:1122490F0
# chr4:163850:163949:-:163911#/1	ATTATTGTAGATTATGTAATTAATTGGTTTTTTAGGTTCTGTTTTTTTTTTGAGTTAGTT	1:0:0:0	RC_chr4:1415973F0
# chr4:834155:834254:+:834178:834192#/1	TTTTATTAGAATAATTTGGAAATCTAGTAGTTAAGAGCTTAAGAGGGTGTAAGGTATGGT	1:0:0:0	chr4:834155F0
sub reads_match{
    my ($length, $bowtie_file, $fr) = @_;
    open my $fh, '<', $bowtie_file;
    while (defined(my $line = <$fh>)){
        chomp $line;
        my ($read_id, $read_seq, $mismatch_code, $target_code) = split /\t/, $line;
        next if $mismatch_code eq 'NM';
        if ($target_code =~ /^(RC_)?(\w+?):(\d+)(F|R)/){
            my ($rc, $seq, $pos) = ($1, $2, $3);
            my $original = $fr->get($seq, $pos, $pos + $length - 1, rc => $rc, coord => ($rc?'r':'f'));
            if (numdiff(c2t($original), c2t($read_seq)) > 3){
                return 0;
            }
        }
        else {
            die "invalid target code $target_code";
        }
    }
    close $fh;
    return 1;
}

sub junk_count{
    my $bt = shift;
    my $NMcount = 0;
    open my $fh, '<', $bowtieout;
    while (defined(my $line = <$fh>)){
        chomp $line;
        my @split = split /\t/, $line;
        if ($split[2] eq 'NM'){
            $NMcount++;
        }
    }
    close $fh;
    return $NMcount;
}

