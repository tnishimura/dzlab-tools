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

my $result = GetOptions ( 
    "reference|r=s" => \my $reference,
    "prefix|p=s"    => \(my $prefix),
    "max=i"         => \(my $max = 50),
);
pod2usage(-verbose => 2, -noperldoc => 1) 
if (!$result || ! $reference || ! -f $reference || ! $prefix );  

my $output = $prefix . ".txt";
my $svg    = $prefix . ".svg";
my $gp     = $prefix . ".gp";

say STDERR "reading in $reference";
my $fr = FastaReader->new(file => $reference);
my $methsites = $fr->num_methylation_sites();
say STDERR "$reference has $methsites methsites";

# methyl-coverage.pl -r reference.fas -- --LANE1 LANE1/single-c/*.gff --LANE2 LANE2/single-c/*.gff
# => {
#   LANE1 => LANE1/single-c/*.gff,
#   LANE2 => LANE2/single-c/*.gff,
# }
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
say STDERR Dumper \%files;

#######################################################################

open my $output_fh, '>', $output;

for my $title (keys %files) {
    # body...
    my $line_count = 0;
    my %coverage_count;

    for my $gff_file (@{$files{$title}}){
        next if ($gff_file =~ /chrm|chrc|chrpt|chrunknown/i);

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
    say $output_fh "# $title";
    for my $c (sort { $a <=> $b } keys %coverage_count) {
        say $output_fh "$c\t$coverage_count{$c}";
    }
    # two blank lines for gnuplot index
    print $output_fh "\n\n";
}
close $output_fh;

#######################################################################
# gnuplot

# read .txt file to create gnuplot file and svg file

open my $gp_fh, '>', $gp;
{
    my $i = 0;
    my @titles = keys %files;

    my $cmd = qq{
    set log y
    set terminal svg size 1500, 750
    set output "$svg"
    plot [0:50] } .  join ", ", map { qq{"$output" i $_ u 1:2 w lp t "$titles[$_]" } } (0 .. $#titles);

    say STDERR $cmd;

    say $gp_fh $cmd;
    system("gnuplot $gp");
}

close $gp_fh;

=head1 NAME

 methyl-coverage.pl -r reference.fas -p prefix -- --LANE1 LANE1/single-c/*.gff --LANE2 LANE2/single-c/*.gff
 
=cut

__END__
set terminal svg size 1500, 750
set output "test.svg"
plot [0:50] "coverage.txt" i 0 u 1:($2 / 162152310) w lp t "10%", \
            "coverage.txt" i 1 u 1:($2 / 162152310) w lp t "20%", \
            "coverage.txt" i 2 u 1:($2 / 162152310) w lp t "40%", \
            "coverage.txt" i 3 u 1:($2 / 162152310) w lp t "60%", \
            "coverage.txt" i 4 u 1:($2 / 162152310) w lp t "80%"

