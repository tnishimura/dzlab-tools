#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use FastaReader;

use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) if $opt_help;

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}

say $opt_sequence;
say $opt_range{'start'};
say $opt_range{'end'};

my $r = FastaReader->new(file => $opt_fasta);

say Dumper $r->length;
say $r->get($opt_sequence, $opt_range{start}, $opt_range{end}, coord => 'f', base => 1);

=head1 NAME
 
fasta.pl
 
=head1 SYNOPSIS

 fasta.pl

=head1 REQUIRED ARGUMENTS

=over

=item  -s <seq> | --sequence <seq>

=item  -r <start> <end> | --range <start> <end>

=for Euclid
    start.type:        int, start >= 1 
    start.type.error:  <fasta> must be > 1
    end.type:        int, end >= 1 
    end.type.error:  <fasta> must be > 1

=item  <fasta> 

=for Euclid
    fasta.type:        readable

=back

=head1 OPTIONS

=over


=item -o <file> | --output <file>

=for Euclid
    file.default:     '-'

=item -h | --help

=back

=cut

