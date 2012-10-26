# this deprecates ParserRole.pm, b/c if you provide a BUILD in a role and you
# consume it in a class that has its own BUILD, the role's BUILD is not run
# b/c Moose can't determine the order to run them in

# this is a Moose base class to encapsulate parser boilerplate code.
# The consumer of the role should provide a next() method which 
# reads an entry from $self->filehandle().
# The end user should do something like:
# my $p = Parser->new(file => $file_or_filehandle)
# while (defined(my $entry = $p->next())){
#     body...
# }

package Parser;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Moose;
use Carp;
use autodie;    
use DZUtil qw/open_filename_or_handle/;

has filename_or_handle => (
    is => 'ro',
    required => 1,
    init_arg => 'file',
);

has filehandle => (
    is => 'rw',
    init_arg => undef,
);

has file => (
    is => 'rw',
    init_arg => undef,
);

# note the :crlf PerlIO layer!

sub BUILD{
    # warn "parser role constructor called";
    my ($self) = @_;
    my ($file, $fh) = open_filename_or_handle($self->filename_or_handle);

    $self->file($file);
    $self->filehandle($fh);
}

sub DEMOLISH{
    my ($self) = @_;
    if (defined $self->file()){
        close $self->filehandle;
    }
}

1;

