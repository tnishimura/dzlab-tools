#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../Tree-Range/lib";
use GFF::Tree;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use Log::Log4perl qw/:easy/;
use Sam::Parser;

Log::Log4perl->easy_init({ level => $DEBUG, layout => '%d{HH:mm:ss} %.1p > %m%n' });
my $logger = get_logger();

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || !$opt_gff || !$opt_input || !$opt_tag;

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}

$logger->info("reading $opt_gff into GFF::Tree...");
my $gt = GFF::Tree->new(file => $opt_gff,normalize => 1, lenient => 1);
$logger->info("done");

my %alignments = do{
    my $pp = GFF::Parser->new(file => $opt_gff);
    my %tmp;
    while (defined(my $gff = $pp->next())){
        if (defined(my $locus = $gff->get_column($opt_tag))){
            $tmp{$locus}=[uc $gff->sequence,$gff->start,$gff->end,$gff->strand,0];
        }
    }
    %tmp;
};

my $counter = 0;

my $sam_reader = Sam::Parser->new(file => $opt_input, skip_unmapped => 1);
while (defined(my $sam = $sam_reader->next())){
    my $strand = $sam->is_reverse() ? '-' : '+';
    my ($seq, $start, $end) = (uc $sam->seqid, $sam->leftmost, $sam->rightmost);

    my @results = $gt->search_overlap($seq,$start,$end);

    for my $result (@results){
        if (defined(my $locus = $result->{item}->get_column($opt_tag))){
            $alignments{$locus}[4]++;
        }
        else{
            die sprintf("%s does not have an $opt_tag field?", $result->{item}->to_string);
        }
    }
    $counter++;
    if ($counter % 100_000 == 0){
        warn $counter;
    }
}

my $num_million_mapped_reads = $counter / 1_000_000;

for my $id (sort keys %alignments) {
    my ($seq, $start, $end, $strand, $n) = @{$alignments{$id}};
    my $size_in_kilo_bases = $end - $start + 1;
    my $rpkm = $n / ($size_in_kilo_bases * $num_million_mapped_reads);

    say join "\t", $seq, '.', $opt_feature, $start, $end, $rpkm, '.', '.', "n=$n";
}

=head1 NAME

 window_by_annotation.pl -g annotation.gff -k -t 'ID' -o output.gff input.sam

=head1 OPTIONS

=over


=item  -g <file> | --gff <file>

Annotation file.

=for Euclid
    file.type:        readable

=item  <input>

=for Euclid
    input.type:        readable

=item  -t <tag> | --tag <tag>

Locus tag in --gff annotation file. Defaults to 'ID'.

=item -o <file> | --output <file>

=for Euclid
    file.default:     '-'

=item  -f <featurename> | --feature <featurename>

=for Euclid
    featurename.default:     'window'

=item  -k | --no-skip 

=item --help | -h

=back

=cut



