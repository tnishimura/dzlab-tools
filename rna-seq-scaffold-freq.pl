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
if !$opt_scaffold || !$opt_left_reads;

(my $prefix = $opt_scaffold) =~ s/\.\w+$//;
$prefix=basename($prefix);
(my $bowtie = $opt_left_reads) =~ s/\.\w+$/.vs-scaffold_$prefix.bowtie/;
my $output = $opt_output // $bowtie . ".freq";
my $log = "$bowtie.log";

say $opt_left_reads;
say $bowtie;


if (defined $opt_right_reads){
    launch("bowtie $opt_scaffold -1 $opt_left_reads -2 $opt_right_reads ?? -v 2 -B 1 --best",
        expected => $bowtie,
        also => $log,
        dryrun => $opt_dry,
    );
}
else {
    launch("bowtie $opt_scaffold $opt_left_reads ?? -v 2 -B 1 --best",
        expected => $bowtie,
        also => $log,
        dryrun => $opt_dry,
    );
}

launch(qq{perl -S parse_bowtie.pl -i "$opt_regex"  -r $opt_scaffold -o ?? -f $bowtie},
    expected => $output,
    dryrun => $opt_dry,
);

=head1 OPTIONS

=over

=item  -s <fasta> | --scaffold <fasta>

=for Euclid
    fasta.type:        readable

=item  -i <regex> | --regex <regex>

=for Euclid
    regex.default:     '([^\s]+)'

=item  -o <file> | --output <file>

=item  -n | --dry 

=item -1 <input> | --left-reads <input>

=for Euclid
    input.type:        readable

=item -2 <input> | --right-reads <input>

=for Euclid
    input.type:        readable

=back

=cut
