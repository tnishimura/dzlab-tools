#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use lib "$FindBin::Bin/Tree-Range/lib";
use FastaReader;
use GFF::Parser;
use GFF::Tree;
use List::Util qw/sum min max/;
use Pod::Usage;
use Getopt::Long;
END {close STDOUT}
$| = 1;

my $opt_window_size = 50;
my $opt_no_skip = 1;
my $opt_num_headers = 5;
my $opt_feature_name = 'part';

my $result = GetOptions (
    "feature|f=s" => \$opt_feature_name,
);
my ($annotation, $reference) = @ARGV;

pod2usage(-verbose => 1) if (!$result);  
if (!$result || ! $annotation || !$reference){
    say "\nusage: te-explode.pl -f feature-name annotation.gff reference.fasta";
    say <<'END';

Description:
Take an annotation, a reference, and an option feature name. For each 50 bp
window, assign score according to its location relative to annotation element,
such that windows aligning to 0-50 of a TE are given .2, 51-100 are given .4,
etc, up to 1.0:

0    50   100  150  200  250  ...
|----|----|----|----|----|----|---...
  .2   .4   .6   .8   1.0  1.0  1.0 ...

If window doesn't properly align to a 50 bp grid relative to the ends (ie, if
it aligns 15-65 of a annotation element), it is given a piecewise score, where
each bp aligning to 1-50 are given a score of .004 ( = .2/50), 51-100 are given
.002 ( = .4/50), and so on. This is not ideal, since it's discontinuous.

END
    exit 1;
}

#######################################################################

my $tree = GFF::Tree->new(file => $annotation, normalize => 1, lenient => 1);
my $fr = FastaReader->new(file => $reference);

for my $seq (sort $fr->sequence_list){
    my $len = $fr->get_length($seq);
    my $start = 1;
    while ($start <= $len){
        my $end = $start + $opt_window_size - 1;
        my @overlaps = $tree->search(uc $seq, $start, $end - 1);
        my $score = scalar @overlaps ?  window_score(\@overlaps, $start, $end) : '.';
        my $attr = scalar @overlaps ?  
        ('ID=' . join ",", map { $_->get_column('ID') } @overlaps)
        : '.';
        if ($score ne '.' || $opt_no_skip){
            say join "\t", $seq, q{.}, $opt_feature_name, $start, $end, $score, qw{. .}, $attr;
        }
        $start += $opt_window_size;
    }
}

sub make_upstream_bp_score_function{
    my $gff = shift;
    my $start = $gff->start;
    return sub{
        my $x = shift;
        given ($x){
            when ($_ < $start){ return 0; }
            when (between($start, $x, $start + 50)){ return .004; }
            when (between($start + 50 + 1, $x, $start + 100)){ return .008; }
            when (between($start + 100 + 1, $x, $start + 150)){ return .012; }
            when (between($start + 150 + 1, $x, $start + 200)){ return .016; }
            when ($start + 201 <= $_) { return .020; }
            default { die "bad arg"; }
        }
    };
}

sub make_downstream_bp_score_function{
    my $gff = shift;
    my $end = $gff->end;
    return sub{
        my $x = shift;
        given ($x){
            when ($_ < $end - 200){ return 0.20; }
            when (between($end - 200, $x, $end - 150)){ return .016; }
            when (between($end - 150 + 1, $x, $end - 100)){ return .012; }
            when (between($end - 100 + 1, $x, $end - 50)){ return .008; }
            when (between($end - 50 + 1, $x, $end - 0)){ return .004; }
            when ($end <= $_) { return .0; }
            default { die "bad arg"; }
        }
    };
}

sub make_bp_score_function{
    my $gff = shift;
    my $u = make_upstream_bp_score_function($gff);
    my $d = make_downstream_bp_score_function($gff);
    return sub {
        my $x = shift;
        return min($u->($x), $d->($x));
    }
}
sub window_score{
    my ($gff_aref, $start, $end) = @_;
    my @accum;
    for my $gff (@$gff_aref) {
        my $f = make_bp_score_function($gff);
        #say Dumper $gff;
        #say Dumper $f;
        #say $start;
        #say $end;
        #say $f->(17050);
        push @accum, (sum map { $f->($_) }  ($start .. $end));
    }
    #die Dumper \@accum;
    return sprintf("%.4f", max(@accum));
}

sub between{
    my ($low, $x, $high) = @_;
    if ($low <= $x && $x <= $high){
        return 1;
    }
    else {
        return 0;
    }
}
