#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use lib "$FindBin::Bin/Tree-Range/lib";
use GFF::Tree;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use Log::Log4perl qw/:easy/;
use Counter;
Log::Log4perl->easy_init({ level => $DEBUG, layout => '%d{HH:mm:ss} %.1p > %m%n' });
my $logger = get_logger();

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || !$opt_gff || !$opt_input || !$opt_tag;

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}

$logger->info("reading $opt_gff into GFF::Tree...");
my $gt = GFF::Tree->new(file => $opt_gff,normalize => 1);
$logger->info("done");

my $p = GFF::Parser->new(file => $opt_input,normalize => 1);

my %methylation = (); # { id => [seq, #c, #t, #n] }

my $counter = Counter->new();

while (defined(my $gff = $p->next())){
    my ($seq, $start, $end, $strand, $c, $t) 
    = ($gff->sequence, $gff->start, $gff->end, $gff->strand, $gff->get_column('c'), $gff->get_column('t'),);
    $strand //= '.';

    my $gffstring = $gff->to_string;

    # if ($c+$t==0){ die "$gffstring"; }
    my @results = $gt->search_overlap($seq,$start,$end);

    for my $result (@results){
        my $overlap = $result->{overlap};
        #say $gff->start . " $overlap ". $overlap/$gff->length();

        if (defined(my $locus = $result->{item}->get_column($opt_tag))){
            if (! exists $methylation{$locus}){
                $methylation{$locus}=[$seq,$start,$end,$strand,0,0,0];
            }
            #my $metharray = $methylation{$locus};
            $methylation{$locus}[4]+=$c;
            $methylation{$locus}[5]+=$t;
            $methylation{$locus}[6]+=1;
        }
        else{
            die sprintf("%s does not have an $opt_tag field?", $result->{item}->to_string);
        }
    }
    $counter->increment();
}

for my $id (sort keys %methylation) {
    my ($seq, $start, $end, $strand, $c, $t, $n) = @{$methylation{$id}};
    #if ($c+$t==0){ die "if \$c+\$t is 0 then why is there an entry?"; }
    my $score = $c+$t == 0 ? 0 : sprintf("%.4f", $c/($c+$t));
    say join("\t", $seq, 'win', $opt_feature, $start, $end, $score, $strand, '.',
        "ID=$id;c=$c;t=$t;n=$n");
}



=head1 NAME

overlaps_gff.pl 

=head1 SYNOPSIS

Usage examples:

 overlaps_gff.pl -g annotation.gff -k -t 'ID' -o output.gff input.gff

=head1 REQUIRED ARGUMENTS

=over

=back

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

=item --help | -h

=back

=cut



