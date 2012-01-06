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
if !$opt_reference || !$opt_left_reads;

(my $prefix = $opt_reference) =~ s/\.\w+$//;
$prefix = basename($prefix);
(my $bowtie = $opt_left_reads) =~ s/\.\w+$/.vs-genome_$prefix.bowtie/;
my $gff = $bowtie . ".gff";
my $w50 = $bowtie . ".w50.gff";
my $log = $bowtie . ".w50.gff.log";

my $fasta_flag = $opt_fasta ? ' -f ' : '';
my $thread_flag = $opt_threads ? "-p $opt_threads" : "";

if (defined $opt_right_reads){
    launch("bowtie $opt_reference -1 $opt_left_reads -2 $opt_right_reads ?? -v 2 -B 1 --best $fasta_flag $thread_flag",
        expected => $bowtie,
        also => $log,
        dryrun => $opt_dry,
    );
}
else {
    launch("bowtie $opt_reference $opt_left_reads ?? -v 2 -B 1 --best $fasta_flag $thread_flag",
        expected => $bowtie,
        also => $log,
        dryrun => $opt_dry,
    );
}
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

=item -1 <input> | --left-reads <input>

=for Euclid
    input.type:        readable

=item -2 <input> | --right-reads <input>

=for Euclid
    input.type:        readable

=item  -f | --fasta 

Reads are FASTA, not FASTQ

=item -p <threads> | --threads <threads>

=back

=cut
