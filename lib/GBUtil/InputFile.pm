package GBUtil::InputFile;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Moose::Role;
use Carp;
use autodie;    

has [qw/file staging_dir source/] => (
    is => 'ro',
    required => 1,
);

has staging_file => (
    is => 'rw',
);

sub title{
    my $self = shift;
    sprintf("%s [%s]", $self->file, $self->source);
}

requires 'convert';
requires 'upload_files';

1;
