#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use List::Util qw/sum/;
use Cwd qw/abs_path/;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
unless ($opt_label xor $opt_output) && $opt_ecotype_a && $opt_ecotype_b && 
$opt_bowtie_mismatches && $opt_eland_filtered_b && $opt_eland_filtered_a; 

use Log::Log4perl qw/:easy/;
Log::Log4perl->easy_init( { 
    level    => $INFO,
    #file     => ">run.log",
    layout   => '%d{HH:mm:ss} %p> (%L) %M - %m%n',
} );
my $logger = get_logger();

#######################################################################
# Count and calc ratios and stuff

my %counts = ($opt_ecotype_a => {}, $opt_ecotype_b => {});
my @chromosomes = qw/ chr3 chr1 chr4 chrc chr2 chrm chr5/;
my @core_chromosomes = qw/ chr3 chr1 chr4 chr5 chr2/;
for my $c (@chromosomes) {
    for my $mm (0..$opt_bowtie_mismatches, 'total') {
        $counts{$opt_ecotype_a}{$c}{$mm} = 0;
        $counts{$opt_ecotype_b}{$c}{$mm} = 0;
    }
}
$logger->debug(Dumper \%counts);

for my $eco ([$opt_ecotype_a, $opt_eland_filtered_a], [$opt_ecotype_b, $opt_eland_filtered_b]){
    my ($ecotype, $file) = @$eco;
    open my $fh, '<', $file or die "can't open $file";
    while (defined(my $line = <$fh>)){
        my @split = split /\t/, $line;
        $split[3] =~ /(?:RC_)?(chr.*?):.*?(\d)$/i;
        my $c = lc $1;
        $counts{$ecotype}{$c}{$2}++;
        $counts{$ecotype}{$c}{total}++;
    }
    close $fh;
}

$logger->debug(Dumper \%counts);

my $ratio_log = $opt_output || "$opt_label.$opt_ecotype_a-vs-$opt_ecotype_b.log";
open my $ratiofh, '>', $ratio_log;

say $ratiofh "$opt_ecotype_a = " . abs_path($opt_eland_filtered_a);
say $ratiofh "$opt_ecotype_b = " . abs_path($opt_eland_filtered_b);

say $ratiofh "\n$opt_ecotype_a / $opt_ecotype_b";

for my $mm (0 .. $opt_bowtie_mismatches, 'total') {
    say $ratiofh "$mm mismatches";
    for my $c (sort @chromosomes) {
        say $ratiofh "$c: \t $counts{$opt_ecotype_a}{$c}{$mm} / $counts{$opt_ecotype_b}{$c}{$mm} = "  . 
        ( $counts{$opt_ecotype_b}{$c}{$mm} > 0 ?
            ($counts{$opt_ecotype_a}{$c}{$mm} / $counts{$opt_ecotype_b}{$c}{$mm})
            : "inf" );
    }
    my $total_a = sum map { $counts{$opt_ecotype_a}{$_}{$mm} } @core_chromosomes;
    my $total_b = sum map { $counts{$opt_ecotype_b}{$_}{$mm} } @core_chromosomes;
    say $ratiofh "total (no c/m): \t $total_a / $total_b = " . ( $total_b > 0 ?  $total_a / $total_b  : 'inf');
}


=head1 NAME

split_ratio.pl - generate stats about files generated by split_on_mismatches.pl

=head1 SYNOPSIS

Usage examples:

 split_ratio.pl -l s3 -ea Col -eb Ler -a split-a.eland -b split-b.eland -o ratio.log

=head1 OPTIONS

=over

=item  -l <label> | --label <label>

Label.

=for Euclid
    label.default:     ''

=item  -o <name> | --output <name>

output file name

=for Euclid
    name.default:     ''

=item  -ea <eco> | --ecotype-a <eco>

Ecotype A label.

=item  -eb <eco> | --ecotype-b <eco>

Ecotype B label.

=item  -a <eland> | --eland-filtered-a <eland>

Post-split_on_mismatch.pl eland A

=for Euclid
    eland.type:        readable

=item  -b <eland> | --eland-filtered-b <eland>

Post-split_on_mismatch.pl eland B

=for Euclid
    eland.type:        readable

=item  -m <num> | --bowtie-mismatches <num>


Number of mismatches to allow in bowtie

=for Euclid
    num.default:     2
    num.type:        int, num >= 0 && num <= 3

=item --help | -h

=back

=cut

