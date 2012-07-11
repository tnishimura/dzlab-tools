#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use File::Basename qw/basename dirname/;
use File::Spec::Functions qw/rel2abs canonpath catdir catfile updir/;
use YAML qw/Load Dump LoadFile DumpFile/;
use Pod::Usage;
use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use FastaReader;

END {close STDOUT}
$| = 1;

my $result = GetOptions ( "reference|r=s" => \my $reference,);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || ! $reference || ! -f $reference);  

say STDERR "reading in $reference";
my $fr = FastaReader->new(file => $reference);
my $methsites = $fr->num_methylation_sites();
say STDERR "$reference has $methsites methsites";

use Tie::IxHash;
tie my %files, 'Tie::IxHash'; 
{
    my $current_key;
    for (@ARGV) {
        if (s/^--//){
            $current_key = $_;
        }
        elsif (defined $current_key){
            push @{$files{$current_key}}, $_;
        }
        else{
            pod2usage(-verbose => 2, -noperldoc => 1);
        }
    }
}

for my $title (keys %files) {
    # body...
    my $line_count = 0;
    my %coverage_count;

    for my $gff_file (@{$files{$title}}){
        next if ($gff_file =~ /chrm|chrc/i);

        say STDERR "processing $gff_file";

        #my ($wcl) = `wc -l $gff_file` =~ /^(\d+)/;
        #$line_count += $wcl;

        my ($context) = basename($gff_file) =~ /\b(CG|CHG|CHH)\b/;
        die "can't parse $gff_file for context" if ! defined $context;

        open my $input_fh, '<', $gff_file;
        #my $p = GFF::Parser->new(file => $gff_file);
        #while (defined(my $gff = $p->next())){
            #my ($c, $t) = ($gff->get_column('c'), $gff->get_column('t'));
        while (defined(my $line = <$input_fh>)){
            my ($c, $t) = $line =~ /c=(\d+);t=(\d+)/;

            if (defined $c and defined $t){
                if ($c + $t == 0){
                    die "\$c + \$t == 0 in $gff_file";
                }
                ++$coverage_count{$c + $t};
                ++$line_count;
            }
            else{
                die "no c or t? then why is it in single-c file?";
            }
            say STDERR $line_count if $line_count % 250000 == 0;
        }
        close $input_fh;
    }

    say STDERR "done with gff files, adding uncounted methyl sites";

    $coverage_count{0} += $methsites - $line_count;
    say "# $title";
    for my $c (sort { $a <=> $b } keys %coverage_count) {
        say "$c\t$coverage_count{$c}";
    }
    # two blank lines for gnuplot index
    print "\n\n";
}

=head1 NAME

coverage.pl -r reference -- --Meow Meow/single-c/*.gff --Bark Bark/single-c/*.gff

=cut
