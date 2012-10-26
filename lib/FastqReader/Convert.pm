package FastqReader::Convert;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use FastqReader;
use DZUtil qw/c2t g2a reverse_complement open_filename_or_handle/;

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(
    fastq_to_fastq_c2t fastq_to_fastq_g2a fastq_to_fastq_c2t_rc fastq_to_fastq_g2a_rc
    fastq_to_fasta_c2t fastq_to_fasta_g2a fastq_to_fasta_c2t_rc fastq_to_fasta_g2a_rc
    fastq_to_fasta
);

sub _convert{
    my ($from_fasta, $to_fasta, $reverse_comp, $methyl_pattern, $input_file_or_filehandle, $output_file) = @_;
    open my $outfh, '>', $output_file;;

    if (defined($methyl_pattern) and ($methyl_pattern ne 'c2t' and $methyl_pattern ne 'g2a')){
        croak "bs_fastq expects methyl_pattern to be c2t or g2a";
    }
    if (! defined $methyl_pattern && ! $reverse_comp){
        croak "why are you converting from fastq to fastq without any transformation?" if (!$from_fasta && !$to_fasta);
        croak "why are you converting from fasta to fasta without any transformation?" if ($from_fasta && $to_fasta);
    }

    if (! defined $methyl_pattern and ! $reverse_comp and ! $to_fasta){
        croak "why are you converting from fastq to fastq without any transformation?";
    }

    my $fqr = FastqReader->new(file => $input_file_or_filehandle, fasta => $from_fasta);
    while (defined(my $r = $fqr->next())){
        my ($readid, $seq, $qual) = @$r;
        $seq = ! defined $methyl_pattern && ! $reverse_comp ? $seq 
             : ! defined $methyl_pattern &&   $reverse_comp ? reverse_complement($seq) 
             : $methyl_pattern eq 'c2t'  && ! $reverse_comp ? c2t($seq)
             : $methyl_pattern eq 'g2a'  && ! $reverse_comp ? g2a($seq)
             : $methyl_pattern eq 'c2t'  &&   $reverse_comp ? c2t(reverse_complement($seq))
             : $methyl_pattern eq 'g2a'  &&   $reverse_comp ? g2a(reverse_complement($seq))
             : croak "BUG";
        if ($to_fasta){
            say $outfh $readid;
            say $outfh $seq;
        }
        else{
            say $outfh $readid;
            say $outfh $seq;
            say $outfh "+";
            say $outfh $qual;
        }
    }

    close $outfh
}

sub fastq_to_fastq_c2t    { _convert(0, 0, 0, 'c2t', @_); }
sub fastq_to_fastq_g2a    { _convert(0, 0, 0, 'g2a', @_); }
sub fastq_to_fastq_c2t_rc { _convert(0, 0, 1, 'c2t', @_); }
sub fastq_to_fastq_g2a_rc { _convert(0, 0, 1, 'g2a', @_); }
sub fastq_to_fasta_c2t    { _convert(0, 1, 0, 'c2t', @_); }
sub fastq_to_fasta_g2a    { _convert(0, 1, 0, 'g2a', @_); }
sub fastq_to_fasta_c2t_rc { _convert(0, 1, 1, 'c2t', @_); }
sub fastq_to_fasta_g2a_rc { _convert(0, 1, 1, 'g2a', @_); }
sub fastq_to_fasta        { _convert(0, 1, 0, undef, @_); }

sub fasta_to_fasta_c2t    { _convert(1, 1, 0, 'c2t', @_); }
sub fasta_to_fasta_g2a    { _convert(1, 1, 0, 'g2a', @_); }
sub fasta_to_fasta_c2t_rc { _convert(1, 1, 1, 'c2t', @_); }
sub fasta_to_fasta_g2a_rc { _convert(1, 1, 1, 'g2a', @_); }

1;
