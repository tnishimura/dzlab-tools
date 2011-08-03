#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use eow;

END {close STDOUT}
$| = 1;

use Log::Log4perl qw/:easy/;
Log::Log4perl->easy_init({ level => $DEBUG, layout => '%d{HH:mm:ss} %.1p > %m%n' });
my $logger = get_logger();

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "singlec=s" => \(my $singlec),
    "original=s" => \(my $original),
);
pod2usage(-verbose => 1) if (!$result);  

my $counter = counter(20000);
my %original = (); # { "seq~pos" => base }

$logger->info("indexing original");
{
    my $p = GFF::Parser->new(file => $original);
    while (defined(my $gff = $p->next())){
        $original{make_key($gff)} = $gff->source;
        $counter->();
    }
}
$counter->('reset');

$logger->info("reading single-c");
{
    my $total = 0;
    my $notok = 0;
    my %stats = ( cok => 0, tok => 0, cnotok => 0, tnotok => 0, both => 0, neither => 0, nonexistant => 0, total => 0);

    my $p = GFF::Parser->new(file => $singlec);
    while (defined(my $gff = $p->next())){
        $counter->();

        my $key           = make_key($gff);
        my $exists        = exists $original{$key};
        my $original_base = $original{$key};
        ++$stats{total};

        if (! $exists){
            $logger->logwarn("non-existant meth found: " . $gff->to_string);
            ++$stats{nonexistant};
            next;
        }

        given (methylation($gff)){
            when ('C'){
                if ($original_base eq 'C'){ ++$stats{cok}; }
                else{ ++$stats{cnotok}; }
            }
            when ('T'){
                if ($original_base eq 'T'){ ++$stats{tok}; }
                else{ ++$stats{tnotok}; }
            }
            when ('neither'){
                ++$stats{neither};

            }
            when ('both'){
                ++$stats{both};

            }
            default {
                die "impossible";
            }
        }
    }
    say "cok: $stats{cok}";
    say "tok: $stats{tok}";
    say "cnotok: $stats{cnotok}";
    say "tnotok: $stats{tnotok}";
    say "both: $stats{both}";
    say "neither: $stats{neither}";
    say "nonexistant: $stats{nonexistant}";
    say "total: $stats{total}";
}

sub make_key{
    my $gff =shift;
    return $gff->sequence . "~" . $gff->start;
}

sub methylation{
    my $gff =shift;
    my $c = $gff->get_column('c');
    my $t = $gff->get_column('t');

    if ($c > 0 && $t ==0){
        return 'C';
    }
    elsif ($c == 0 && $t > 0){
        return 'T';
    }
    elsif ($c == 0 && $t == 0){
        return 'neither';
    }
    else {
        return "both";
    }
}
