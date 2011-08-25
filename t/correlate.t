package correlate::Test;
use base qw(Test::Class);
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Test::More 'no_plan';
use FindBin;
use lib "$FindBin::Bin/lib";
use Launch;
use FastaReader;
use DZUtil qw/c2t fastq_convert_read_header/;

our $left = 't/data/bs-sequel-test.fastq_1-60.eland3.post';
our $right = 't/data/bs-sequel-test.fastq_61-100.eland3.post';
our $reads = 't/data/bs-sequel-test.fastq';
our $ref = 't_intermediate/TAIR_mini.fas';
our $singleout = 't_intermediate/correlateSingle.test';
our $pairedout = 't_intermediate/correlatePaired.test';

sub setup : Test(setup){
    my $self = shift;
    $self->{fr} = FastaReader->new(file => $ref,slurp => 1);
    $self->{reads} = slurp_reads($reads);
}

sub test_correlatePairedEnds : Tests{
    my $self = shift;

    ok(launch("perl -S correlatePairedEnds.pl -l $left -r $right -ref $ref -o $pairedout -t 0 -d 0 -s 100 -2 0 -1 1 -m 10 -a 1"),
        "it runs");
    ok(targets_match_genome($pairedout, $self->{fr}, $self->{reads}), 'targets match genome');
    ok(reads_match_fastq($pairedout, $self->{fr}, $self->{reads}), 'reads match fastq');
}
sub test_correlateSingleEnds : Tests{
    my $self = shift;

    ok(launch("perl -S correlateSingleEnds.pl -e $left --reference $ref -o $singleout -m 10 --random-assign 1"),
        "it runs");
    ok(targets_match_genome($singleout, $self->{fr}, $self->{reads}), 'targets match genome');
    ok(reads_match_fastq($singleout, $self->{fr}, $self->{reads}), 'reads match fastq');
}

correlate::Test->runtests;

sub parse_correlate_line{
    my ($line, $fr) = @_;
    chomp $line;
    my ($seq, $read_field, $target_field,$start,$end,$strand) = (split /\t/, $line)[0,2,8,3,4,6];
    return '.' if $seq eq '.';

    my $target = $target_field =~ /target=([A-Z]+)/ ? $1 : die "target malformed";
    my ($read_id, $read) = $read_field =~ /^(.*):([A-Z]+)$/ ? ($1,$2) : die "read malformed";

    my $original = $fr->get($seq, $start-2, $end+2, coord => $strand eq '-' ? 'r' : 'f');

    return ($seq, $start, $end, $strand, $read_id, $read, $target, $original);
}

sub targets_match_genome{
    my ($file, $fr, $reads_href) = @_;

    open my $fh, '<', $file;
    while (defined(my $line = <$fh>)){
        chomp $line;
        my ($seq, $start, $end, $strand, $read_id, $read, $target, $original) = 
        parse_correlate_line($line, $fr);
        next if $seq eq '.';

        if ($original ne $target){
            say STDERR "target mismatch: ($.)\n$line\n$original\n$target";
            return 0;
        }
    }
    close $fh;
    return 1;
}

sub reads_match_fastq{
    my ($file, $fr, $reads_href) = @_;

    open my $fh, '<', $file;
    while (defined(my $line = <$fh>)){
        chomp $line;
        my ($seq, $start, $end, $strand, $read_id, $read, $target, $original) = 
        parse_correlate_line($line, $fr);
        next if $seq eq '.';

        if (! exists $reads_href->{$read_id}){
            say STDERR "no such read $read_id in original (line $.):\n$read";
            return 0;
        } 
        else{
            my $part;
            if (length $read == 60){ # left half
                $part = substr $reads_href->{$read_id},0,60;
            } 
            else{
                $part = substr $reads_href->{$read_id},60,40;
            }
            if (c2t($part) ne c2t($read)){
                say STDERR "read mismatch $read_id (line$.)\n$read\n$part";
                return 0
            }
        }
    }
    close $fh;
    return 1;
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

sub slurp_reads{
    my ($reads) = @_;
    my %accum;
    open my $fh, '<', $reads;
    while (! eof $fh){
        my @lines = map {scalar <$fh>} (0..3);
        chomp @lines;
        $accum{fastq_convert_read_header($lines[0])} = $lines[1];
    }
    close $fh;
    return \%accum;
}
1;
