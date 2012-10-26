package FastqReader::CountReads;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use FastqReader;
use DZUtil qw/c2t g2a open_filename_or_handle/;

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(count_reads count_reads_fasta);

sub _count_reads{
    my ($input_file_or_filehandle, $fasta) = @_;
    my $fqr = FastqReader->new(file => $input_file_or_filehandle, fasta => $fasta);

    my $count = 0;
    while (defined(my $quartet = $fqr->next())){
        $count++;
    }

    return $count;
}

sub count_reads { _count_reads(@_, 0); }
sub count_reads_fasta { _count_reads(@_, 1); }

1;
