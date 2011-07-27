#!/usr/bin/env perl
use 5.10.0;
use strict;
use warnings FATAL => "all";
use Data::Dumper;
use feature 'say';
use autodie;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/lib";
use Launch;
use File::Basename;
use Log::Log4perl qw/:easy/;
Log::Log4perl->easy_init({ level => $DEBUG, layout => '%d{HH:mm:ss} %.1p > %m%n' });
my $logger = get_logger();

END {close STDOUT}

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if !$opt_reference || !$opt_input;

(my $prefix = $opt_reference) =~ s/\.\w+$//;
$prefix = basename($prefix);
(my $bowtie = $opt_input) =~ s/\.\w+$/.vs-genome_$prefix.bowtie/;
my $gff = $bowtie . ".gff";
my $w50 = $bowtie . ".w50.gff";

launch("bowtie $opt_reference $opt_input ?? -v 2 -B 1 --best",
    expected => $bowtie,
    dryrun => $opt_dry,
);
launch(sprintf(q{perl -S parse_bowtie.pl -g -o - %s | perl -wlnaF'\t' -e '$F[3]=$F[4]; print join "\t",@F' > ?? }, $bowtie),
    expected => $gff,
    dryrun => $opt_dry,
);
launch("window_by_fixed.pl -m -w 50 -k -r $opt_reference -o ?? $gff",
    expected => $w50,
    dryrun => $opt_dry,
);

=head1 OPTIONS

=over

=item  -r <fasta> | --reference <fasta>

=for Euclid
    fasta.type:        readable

=item  -n | --dry 

=item <input>

=for Euclid
    input.type:        readable

=back

=cut
