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

$opt_output //= $opt_left_reads;

# get scaffold's basename
(my $scaffold_prefix = $opt_scaffold) =~ s/\.\w+$//;
$scaffold_prefix=basename($scaffold_prefix);

# splice
my ($splice_start, $splice_end) = @opt_splice{'start', 'end'};
my $splicearg;
my $splice_note;
if (defined $opt_length && defined $splice_start && defined $splice_end){
    $splicearg = sprintf("-5 %d -3 %d", $splice_start - 1, $opt_length - $splice_end);
    $splice_note = "$splice_start-$splice_end";
}
elsif (! defined $opt_length && ! defined $splice_start && ! defined $splice_end){
    $splicearg = "";
    $splice_note = "full";
}
else{
    say "if you're going to specifiy splice, you need --splice AND --length";
    exit 1;
}

my $bowtie = "$opt_output.vs-scaffold_$scaffold_prefix.$splice_note.bowtie";
my $output = "$bowtie.freq";
my $output_forward = "$bowtie.forward.freq";
my $output_reverse = "$bowtie.reverse.freq";
my $log = "$bowtie.log";
my $threadarg = $opt_threads ? "-p $opt_threads" : "";
my $fastaflag = $opt_fasta ? '-f' : '';


if (defined $opt_right_reads){
    launch("bowtie $fastaflag $opt_scaffold $splicearg $threadarg -1 $opt_left_reads -2 $opt_right_reads ?? -v $opt_mismatches -B 1 --best",
        expected => $bowtie,
        also => $log,
        dryrun => $opt_dry,
    );
}
else {
    launch("bowtie $fastaflag $opt_scaffold $splicearg $threadarg $opt_left_reads ?? -v $opt_mismatches -B 1 --best",
        expected => $bowtie,
        also => $log,
        dryrun => $opt_dry,
    );
}

if ($opt_split_strand){
    launch(qq{perl -S bowtie_count.pl --forward-only -r $opt_scaffold -o ?? $bowtie},
        expected => $output_forward,
        dryrun => $opt_dry,
    );
    launch(qq{perl -S bowtie_count.pl --reverse-only -r $opt_scaffold -o ?? $bowtie},
        expected => $output_reverse,
        dryrun => $opt_dry,
    );
}
else {
    launch(qq{perl -S bowtie_count.pl -r $opt_scaffold -o ?? $bowtie},
        expected => $output,
        dryrun => $opt_dry,
    );
}


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

=item -p <threads> | --threads <threads>

=item  --splice <start> <end> 

=item -l <len> | --length <len>

=item  -m <mm> | --mismatches <mm>

-v to bowtie

=for Euclid
    mm.default:     2

=item  -f | --fasta

Scaffold is fasta, not a fastq.

=item  --split-strand

=back

=cut
