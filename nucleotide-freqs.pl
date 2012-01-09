#!/usr/bin/env perl
# ___UNDOCUMENTED___

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/lib";
use FastaReader;


# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV;

my $kmer = 2;
my @excluded_sequences;
my $alphabet = 'ACGT';
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'kmer|k=i'     => \$kmer,
    'exclude-seq|x=s{,}' => \@excluded_sequences,
    'alphabet|a=s' => \$alphabet,
    'output|o=s' => \$output,
    'verbose|v'  => sub { use diagnostics; },
    'quiet|q'    => sub { no warnings; },
    'help|h'     => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'   => sub { pod2usage ( -verbose => 2 ); }
);

if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

my $reference = $ARGV[0];
my %nucleotide_frequencies;

print STDERR "# Reading in $reference";
my $fr = FastaReader->new(file => $reference, slurp => 0);
print STDERR "# Done";

for my $k (1, $kmer) {
    
    print STDERR "#k:\t$k\n";

    my %frequencies  = ();
    my $total_length = 0;

    for my $chromosome ($fr->sequence_list()) {
        my $chrlen = $fr->get_length($chromosome);

        $total_length += $chrlen;
        
        print STDERR "#chromosome:\t$chromosome\n";
        print STDERR "#length:\t", $chrlen, "\n";

        my %frequency = %{ word_composition ($fr->get($chromosome, undef, undef), $k) };

        for my $word (sort keys %frequency) {

            next if $word =~ m/[^$alphabet]/i;

            $frequencies{$word} += $frequency{$word};

            print STDERR join ("\t",
                               $word,
                               $frequency{$word},
                               $frequency{$word} / $chrlen, 
                           ), "\n";
        }
    }

    print STDERR "#all:\n";
    print STDERR "#length:$total_length\n";

    print join ("\t",
                '#word',
                '#count',
                '#size',
                '#observed',
                '#expected',
                '#obs/expect',
                '#independent',
            ), "\n"
            if $k > 1;

    for my $word (sort keys %frequencies) {

        next if $word =~ m/[^$alphabet]/i;

        if ($k == 1) {
            $nucleotide_frequencies{$word} = $frequencies{$word} / $total_length;
            print STDERR join ("\t",
                        $word,
                        $nucleotide_frequencies{$word},
                    ), "\n";
        }
        else {
            my $observed = $frequencies{$word} / $total_length;
            my $expected = 1;
            map {$expected *= $nucleotide_frequencies{$_}} (split //, $word);

            print join ("\t",
                        $word,
                        $frequencies{$word},
                        $total_length,
                        $observed,
                        $expected,
                        $observed / $expected,
                        map {$nucleotide_frequencies{$_} * $total_length} (split //, $word),
                    ), "\n";
        }
    }
}



sub word_composition {
    my ($sequence, $k, $alphabet) = @_;

    my %frequency = ();

    for (0 .. length ($sequence) - $k) {
     
        my $word = substr $sequence, $_, $k;
        next if $word =~ m/[^ACGT]/i;
        $frequency{$word}++;

    }

    return \%frequency;
}


__END__


=head1 NAME

 name.pl - Short description

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 OPTIONS

 name.pl [OPTION]... [FILE]...

 -o, --output      filename to write results to (defaults to STDOUT)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page

=head1 REVISION

 Version 0.0.1

 $Rev: 249 $:
 $Author: psilva $:
 $Date: 2010-01-11 21:24:34 -0800 (Mon, 11 Jan 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/nucleotide-freqs.pl $:
 $Id: nucleotide-freqs.pl 249 2010-01-12 05:24:34Z psilva $:

=head1 AUTHOR

 Pedro Silva <psilva@nature.berkeley.edu/>
 Zilberman Lab <http://dzlab.pmb.berkeley.edu/>
 Plant and Microbial Biology Department
 College of Natural Resources
 University of California, Berkeley

=head1 COPYRIGHT

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
