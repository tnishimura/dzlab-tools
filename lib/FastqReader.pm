package FastqReader;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Moose;
use Carp;
use autodie;    
use List::MoreUtils qw/ any /;
use DZUtil qw/c2t g2a open_filename_or_handle/;

extends 'Parser';

has 'fasta' => (
    is => 'rw',
    default => 0,
);

sub eof{
    my $self = shift;
    return eof $self->filehandle;
}

# returns [$readid, $sequence, [$quality]]
# no third line b/c it's just $readid with s/^@/\+/
sub next{
    my $self = shift;
    my $fh = $self->filehandle;

    defined(my $readid = scalar readline $fh) or return;
    defined(my $sequence = scalar readline $fh) or return;
    chomp($readid, $sequence);

    if (! $self->fasta){
        defined(my $readid_again = scalar readline $fh) or return;
        defined(my $quality = scalar readline $fh) or return;
        chomp($readid_again, $quality);

        # $readid =~ s/^@//;
        # $readid_again =~ s/^\+//;
        return [$readid, $sequence, $quality];
    }
    else{
        # $readid =~ s/^@/>/;
        return [$readid, $sequence];
    }
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
