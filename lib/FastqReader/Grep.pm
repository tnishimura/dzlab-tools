package FastqReader::Grep;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(fastq_grep);

sub fastq_grep{
    my $pattern = shift;
    croak "need a pattern" if ! $pattern;
    my $infile = shift;
    my $outfile = shift;

    my $cmd = "grep --no-group-separator -B1 -A2";

    if ($infile && $outfile){
        system("$cmd $pattern $infile > $outfile");
    }
    elsif ($infile && ! $outfile){
        system("$cmd $pattern $infile");
    }
    elsif (! $infile && $outfile){
        system("$cmd $pattern > $outfile");
    }
    else{
        system("$cmd $pattern");
    }
}

1;
