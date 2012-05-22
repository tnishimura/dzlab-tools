#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use Launch;

END {close STDOUT}
$| = 1;

use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || ! $opt_reference || ! $opt_bowtie;

my $gff = $opt_bowtie . ".gff";
my $w50 = $opt_bowtie . ".w50.gff";

launch(sprintf(q{perl -S parse_bowtie.pl -g -o - %s | perl -wlnaF'\t' -e '$F[3]=$F[4]; print join "\t",@F' > ?? }, $opt_bowtie),
    expected => $gff,
    dryrun => $opt_dry,
);
launch("window_by_fixed.pl -m -w 50 -k -r $opt_reference -o ?? $gff",
    expected => $w50,
    dryrun => $opt_dry,
);

=head1 SYNOPSIS

Usage examples:

 window_bowtie.pl 

=head1 OPTIONS

=over

=item -r <fasta> | --reference <fasta>

=for Euclid
    fasta.type:        readable

=item  -w <bases> | --window-size <bases>

=for Euclid
    bases.default:     50

=item <bowtie>

=for Euclid
    bowtie.type:        readable

=item --help | -h

=item  -n | --dry <varname>

=back

=cut
