#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Scalar::Util qw/looks_like_number/;
use FindBin;
use lib "$FindBin::Bin/../lib";
# use Sam::Parser;
use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    #"fix-rc|rc"       => \(my $fixrc),
    "output|o=s" => \(my $output = '-'),
);
if (! $result || (! @ARGV && -t STDIN)){
    pod2usage(-verbose => 2, -noperldoc => 1);
}
my $outfh = $output eq '-' ? \* STDOUT : IO::File->new($output, 'w');

# too slow
# my $parser = Sam::Parser->new(
#     file          => \*ARGV,
#     skip_unmapped => 0,
#     # convert_rc    => ($fixrc // 0),
# );

#######################################################################
# print header lines, but omit irrelevant @SQ lines when filtering by seqid

# my %lengths = %{$parser->length()};

use File::Temp;
my $tmpfh = File::Temp->new( TEMPLATE => 'tempXXXXX', DIR => '.', UNLINK => 1 );
my $tmp_filename = $tmpfh->filename;

warn "tmp filename: $tmp_filename";

use List::Util qw/sum/;
use Scalar::Util qw/looks_like_number/;
# while (defined(my $sam = $parser->next())){ $tmpfh->print($sam->readid() . "\n"); }
while (defined(my $line = <ARGV>)){
    chomp $line;
    next if ($line =~ /^@/);

    my @fields = split /\t/, $line;
    my ($readid, $flags) = @fields;
    my ($mdzstring) = grep { s/^MD:Z://; $_ } @fields[12 .. $#fields];

    die "$ARGV ($.): $line" unless $readid && defined $flags;
    my $mapped = ! ($flags & 0x4);
    my $bases = 0;
    if (defined $mdzstring){
        $bases = sum map { looks_like_number($_) ? $_ : 1 } ($mdzstring =~ /(\d+|\w)/xmg);
    }

    $tmpfh->print("$readid\t$mapped\t$bases\n");
}

$tmpfh->close;

system("sort -S10% -i $tmp_filename -o $tmp_filename");

my $num_reads = 0;
my $num_at_least_one_alignment = 0;
my $num_multiple_alignments = 0;
my $lastid;
my $lastid_already_counted;
my $total_bases;

my $reopened = IO::File->new($tmp_filename);
while (defined(my $line = <$reopened>)){
    chomp $line;
    my ($readid, $mapped, $bases) = split /\t/, $line;
    if (! $mapped){
        $num_reads++;
        next;
    }
    $total_bases += $bases;
    if ($lastid && $readid eq $lastid){
        if ($lastid_already_counted){

        }
        else{
            $num_multiple_alignments++;
            $lastid_already_counted = 1;
        }
    }
    elsif ($readid){
        $num_reads++;
        $num_at_least_one_alignment++;
        $lastid_already_counted = 0;
    }
    $lastid = $readid;
}
$reopened->close;

my $once_percentage = sprintf("%0.2f", 100 * $num_at_least_one_alignment / $num_reads);
my $multiple_percentage = sprintf("%0.2f", 100 * $num_multiple_alignments / $num_reads);
$outfh->print(<<END);
number of reads total: $num_reads
number of reads aligning at least once: $num_at_least_one_alignment ($once_percentage%)
number of reads aligning multiple times: $num_multiple_alignments ($multiple_percentage%)
number of total bases from all aligning reads: $total_bases
END

$outfh->close;

=head1 NAME

Extract alignment statistics (like what bowtie prints to stderr) from sam file.  For this
to be meaningful the sam should have non-mapping reads.

 sam-filter.pl in.sam 

=cut


