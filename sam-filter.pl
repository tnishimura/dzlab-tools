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

#######################################################################
# print header lines, but omit irrelevant @SQ lines when filtering by seqid

for my $hline (@{$parser->header_lines}) {
    # @SQ	SN:chr2	LN:19705359
    # @SQ	SN:chr3	LN:23470805
    if ($sequence_id){
        if ($hline =~ /^\@SQ/){
            if (
                ($fixrc && $hline =~ /\tSN:(?:RC_)?$sequence_id/i)
                ||
                ($fixrc && $hline =~ /\tSN:$sequence_id/i)
            ){
                say $hline;
            }
        }
        else{
            say $hline;
        }
    }
    else{
        say $hline;
    }

}

while (defined(my $sam = $parser->next())){
    next if (defined $min_quality and $sam->mapq < $min_quality);
    my $sam_seqid = $sam->seqid ? $sam->seqid =~ s/^RC_//r : undef;
    next if (defined $sequence_id and defined $sam_seqid and lc($sequence_id) ne lc($sam_seqid));
    say $sam;
}

if ($output){
    close STDOUT;
}

=head1 NAME

 sam-filter.pl [-m|--mapped] [-q|--min-quality 10] [-s|--sequence-id chr1] [--fix-rc|-rc] in.sam > out.sam

=cut

