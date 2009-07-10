#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

my $id_field_name;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'id-field-name|i=s' => \$id_field_name,
    'output|o=s' => \$output,
    'verbose|v'  => sub { use diagnostics; },
    'quiet|q'    => sub { no warnings; },
    'help|h'     => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'   => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result and $id_field_name;


if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

my %genes = ();

while (<>) {
    chomp;
    next if m/^#.*$|^\s*$/;
    my %locus = %{ gff_read ($_) };
    next unless $locus{feature} =~ m/exon/i;

    my ($gene_id) = $locus{attribute} =~ m/.*$id_field_name[\s=]*"?(\w*\d+)"?/;
    croak "Couldn't find the gene id based on the id field name you provided. Either you provided an invalid field name or my expression matching sucks..."
    unless $gene_id;

    my $chr = $locus{seqname};

    $genes{$chr}{$gene_id}{strand} = $locus{strand }  unless exists $genes{$chr}{$gene_id}{strand};
    $genes{$chr}{$gene_id}{start } = $locus{start  }  unless exists $genes{$chr}{$gene_id}{start } and $locus{start} > $genes{$chr}{$gene_id}{start};
    $genes{$chr}{$gene_id}{end   } = $locus{end    }  unless exists $genes{$chr}{$gene_id}{end   } and $locus{start} < $genes{$chr}{$gene_id}{end};
}


for my $chr (sort keys %genes) {

    for my $gene_id (sort { $genes{$chr}{$a}{start} <=>  $genes{$chr}{$b}{start} } keys %{$genes{$chr}}) {

        print join ("\t",
                    $chr,
                    'dz',
                    'gene',
                    $genes{$chr}{$gene_id}{start},
                    $genes{$chr}{$gene_id}{end},
                    q{.},
                    $genes{$chr}{$gene_id}{strand},
                    q{.},
                    "ID=$gene_id",
                    "\n"
                );
    }
}


sub gff_read {
    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, shift);

    $seqname =~ tr/A-Z/a-z/;

    my %rec = (
	'seqname'   => $seqname,
	'source'    => $source,
	'feature'   => $feature,
	'start'     => $start,
	'end'       => $end,
	'score'     => $score,
	'strand'    => $strand,
	'frame'     => $strand,
	'attribute' => $attribute
	);
    return \%rec;
}



__END__


=head1 NAME

 assemble_genes_from_annotation.pl - Build gene models from gff exons annotations

=head1 SYNOPSIS

 assemble_genes_from_annotation.pl -i gene_id_field_name -o gene_annotation.gff no_gene_annotation.gff

=head1 DESCRIPTION

=head1 OPTIONS

 assemble_genes_from_annotation.pl [OPTION]... [FILE]...

 -i, --id-field-name  name immediately preceeding the gene id in the input gff file
 -o, --output         filename to write results to (defaults to STDOUT)
 -v, --verbose        output perl's diagnostic and warning messages
 -q, --quiet          supress perl's diagnostic and warning messages
 -h, --help           print this information
 -m, --manual         print the plain old documentation page

=head1 REVISION

 Version 0.0.1

 $Rev$:
 $Author$:
 $Date$:
 $HeadURL$:
 $Id$:

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