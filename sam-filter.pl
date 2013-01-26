#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Scalar::Util qw/looks_like_number/;
use FindBin;
use lib "$FindBin::Bin/lib";
use Sam::Parser;

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "mapped|m"        => \(my $mapped),
    "min-quality|q=i" => \(my $min_quality),
    "sequence-id|s=s" => \(my $sequence_id),
    "output|o=s"      => \(my $output),
    "fix-rc|rc"       => \(my $fixrc),
);
if (! $result || (! @ARGV && -t STDIN)){
    pod2usage(-verbose => 2, -noperldoc => 1);
}

if ($output){
    open my $fh, '>', $output;
    select $fh;
}

my $parser = Sam::Parser->new(
    file          => \*ARGV,
    skip_unmapped => ($mapped // 0),
    convert_rc    => ($fixrc // 0),
);

say $parser->header;

while (defined(my $sam = $parser->next())){
    next if (defined $min_quality and $sam->mapq < $min_quality);
    next if (defined $sequence_id and $sequence_id ne $sam->seqid);
    say $sam;
}

if ($output){
    close STDOUT;
}

=head1 NAME

sam-only-mapped.pl in.sam > out.sam

=cut

