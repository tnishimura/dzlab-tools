package GBUtil::InputFile::GFF;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use File::Basename qw/basename/;
use File::Spec::Functions qw/catfile/;
use Moose;
use GFF::Parser;

extends 'GBUtil::InputFile';

sub BUILD{
    my $self = shift;
    my %opt = @_;
    $self->(catfile($self->staging_dir, basename($self->file) . ".normalized.gff"));
}

has _features => (
    traits  => ['Array'],
    is      => 'ro',
    isa     => 'ArrayRef[Str]',
    default => sub { [] },
    handles => {
        features => 'elements',
    },
);

# normalize gff. lowercase seqid. fill in source field (default to ".". this matches
# prepare_fasta's metafile's source).
# my $staging_file_name = prepare_gff($input_file_name);
# perl -I$HOME/dzlab-tools/lib/ -MGBUtil -wle 'prepare_gff("foo.bar.gff")'
sub convert{
    my $self = shift;

    my %features;
    open my $outfh, '>', $self->staging_file;
    my $p = GFF::Parser->new(file => $self->file);
    my $source = $self->source;

    while (defined(my $gff = $p->next)){
        $gff->sequence(lc($gff->sequence()) // '.');
        $gff->source($source); 

        if (! exists $features{$gff->feature()}){
            $features{$gff->feature()} = 1;
        }

        say $outfh $gff;
    }
    close $outfh;
    
    $self->_features([sort keys %features]);

    return;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
