# this is a Moose mix-in to encapsulate parser boilerplate code.
# The consumer of the role should provide a next() method which 
# reads an entry from $self->filehandle().
# The end user should do something like:
# my $p = Parser->new(file => $file_or_filehandle)
# while (defined(my $entry = $p->next())){
#     body...
# }

package ParserRole;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Moose::Role;
use Carp;
use autodie;    

requires 'next';

has filename_or_handle => (
    is => 'ro',
    required => 1,
    init_arg => 'file',
);

has filehandle => (
    is => 'rw',
    init_arg => undef,
);

sub BUILD{
    warn "parser role constructor called";
    my ($self) = @_;
    if (ref $self->filename_or_handle eq 'GLOB'){
        $self->filehandle($self->filename_or_handle);
    }
    elsif (!ref $self->filename_or_handle && -f $self->filename_or_handle ){
        open my $fh, '<:crlf', $self->filename_or_handle
            or croak "cannot open " .  $self->filename_or_handle;
        $self->filehandle($fh);
    } elsif (! -f $self->filename_or_handle){
        croak $self->filename_or_handle . " doesn't exist?";
    }
    else {
        croak "file argument to GFF::Parser needs to be file handle or file name" 
        . Dumper $self;
    }
}

sub DEMOLISH{
    my ($self) = @_;
    if (!ref $self->filename_or_handle){
        close $self->filehandle
            or croak "cannot close $self->filename_or_handle";
    }
}

1;
