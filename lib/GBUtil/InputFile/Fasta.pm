package GBUtil::InputFile::Fasta;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use File::Basename qw/basename/;
use File::Spec::Functions qw/catfile/;
use Moose;
use FastaReader;

extends 'GBUtil::InputFile';

sub BUILD{
    my $self = shift;

    $self->staging_file(
        catfile($self->staging_dir, basename($self->file) . ".normalized.fasta")
    );
    $self->meta_file($self->staging_file() . ".meta");
}

has meta_file => (
    is => 'rw',
    init_arg => undef,
);

# normalize fasta header lines, create associated GFF.
# my ($staging_file_name, $meta_file_name) = prepare_fasta($input_file_name);
# perl -I$HOME/dzlab-tools/lib/ -MGBUtil -wle 'prepare_fasta("TAIR_reference.fa")'
sub convert{
    my $self = shift;
    say "+++ " . $self->staging_file;

    # copy input to meta, normalizing headers in the process.
    {
        open my $outfh, '>', $self->staging_file;
        open my $infh, '<:crlf', $self->file;
        while (defined(my $line = <$infh>)){
            chomp $line;
            if ($line =~ /^>(\w+)/){
                say $outfh ">\L$1";
            }
            else{
                say $outfh $line;
            }
        }
        close $infh;
        close $outfh;
    }
    
    # create meta
    my $fr = FastaReader->new(file => $self->staging_file, normalize => 0, slurp => 0);
    my %seqlen = $fr->sequence_lengths;

    open my $metafh, '>', $self->meta_file;
    say $metafh "##gff-version 3\n";
    for my $seqid (sort keys %seqlen) {
        say $metafh join "\t", 
        $seqid, 
        qw/./, # ? 
        'chromosome',
        1,
        $seqlen{$seqid},
        qw/. . ./,
        "ID=$seqid;Name=$seqid",
    }
    close $metafh;

    return;
}

1;
