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
use DZUtil qw/open_filename_or_handle/;

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

    my ($file, $fh) = open_filename_or_handle($filename_or_handle);
    $self->{handle} = $fh;
    $self->{file}   = $file;

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
