# this is a Moose base class to encapsulate parser boilerplate code.
# The consumer of the role should provide a next() method which 
# reads an entry from $self->filehandle().
# The end user should do something like:
# my $p = Parser->new(file => $filename_or_filehandle)
# while (defined(my $entry = $p->next())){
#     body...
# }

package ParserNoMoose;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Hash::Util qw/lock_keys/;
use Params::Validate qw/:all/;
use IO::Scalar;

sub new {
    my $class = shift;
    my $self = bless {}, $class;
    $self->{handle} = undef;
    $self->{file} = undef;

    my %opt = validate(@_, {
            file => {
                callbacks => {
                    'file needs to be filename or a handle' => sub { 
                        my $forh = shift;
                        -f $forh || ref $forh eq 'GLOB' || ref $forh eq 'SCALAR'
                    },
                },
            }, 
            debug => 0,
        });

    my $filename_or_handle = $opt{file};

    # note the :crlf PerlIO layer!
    # handle
    if (ref $filename_or_handle eq 'GLOB'){
        binmode($filename_or_handle, ':crlf');
        $self->{handle} = $filename_or_handle;
    }
    # string
    elsif (ref $filename_or_handle eq 'SCALAR'){ # for debugging
        $self->{handle} = IO::Scalar->new($filename_or_handle);
    }
    # file
    elsif (!ref $filename_or_handle && -f $filename_or_handle ){
        open my $handle, '<:crlf', $filename_or_handle
            or croak "cannot open " .  $filename_or_handle;
        $self->{handle} = $handle;
        $self->{file}   = $filename_or_handle;
    } 
    elsif (! -f $opt{file}){
        croak $filename_or_handle . " #doesn't exist?";
    }
    else {
        croak "file argument to GFF::Parser needs to be file handle or file name" 
        . Dumper $self;
    }

    lock_keys(%$self);
    return $self;
}

sub DESTROY{
    my ($self) = @_;
    if (defined $self->{file}){
        close $self->{handle}
    }
}

1;
