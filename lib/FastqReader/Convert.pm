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
our @EXPORT = qw( fastq_convert);

use Params::Validate qw/:all/;
# types: SCALAR ARRAYREF HASHREF CODEREF GLOB GLOBREF SCALARREF UNDEF OBJECT(blessed) BOOLEAN(UNDEF | SCALAR) HANDLE

sub fastq_convert{
    my %opt = validate(@_, {
            from_fasta     => {
                default => 0,
                optional => 1,
            },
            to_fasta       => {
                default => 1,
                optional => 1,
            },
            rc => 0,
            methyl => {
                callbacks => {
                    'c2t or g2a' => sub { 
                        my $m = shift; 
                        ! defined $m || $m eq 'c2t' || $m eq 'g2a';
                    },
                },
                optional => 1,
            },
            in  => 1,
            out => 1,
        });

    my     ($from_fasta, $to_fasta, $reverse_comp, $methyl_pattern, $input_file_or_filehandle, $output_file ) = 
    @opt{qw/ from_fasta   to_fasta   rc             methyl           in                      out        /};

    my $outfh;
    if (ref $output_file eq 'GLOB'){
        $outfh = $output_file;
    }
    elsif (ref $output_file eq ''){
        open $outfh, '>', $output_file;
    }
    else{
        $outfh = \*STDOUT;
    }

    if (defined($methyl_pattern) and ($methyl_pattern ne 'c2t' and $methyl_pattern ne 'g2a')){
        croak "bs_fastq expects methyl_pattern to be c2t or g2a";
    }
    if ($from_fasta && ! $to_fasta){
        croak "won't convert from a fasta to fastq";
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
            $readid =~ s/^@/>/;
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

    if (-f $output_file){
        close $outfh;
    }
}

1;
