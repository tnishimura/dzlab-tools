package FastqReader;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Moose;
use Carp;
use autodie;    
use List::MoreUtils qw/ any /;


with 'ParserRole';

has 'linesper' => (
    is => 'rw',
    default => 4,
);

# given a list of read ids, return a hash of { id => sequence }
sub get_reads{
    my $self = shift;
    my @query_ids = @_;
    return {} if ! @query_ids;

    my %results;
    while (defined(my $q = $self->next())){
        my ($readid, $sequence) = @$q;
        QUERY:
        for my $query (@query_ids) {
            if ($readid =~ /\Q$query\E/){
                $results{$query} = $sequence;
                last QUERY;
            }
        }
    }
    return \%results;
}

sub next{
    my $self = shift;
    my $fh = $self->filehandle;

    defined(my $readid = scalar readline $fh) or return;
    defined(my $sequence = scalar readline $fh) or return;
    chomp($readid, $sequence);

    if ($self->linesper() == 4){
        defined(my $readid_again = scalar readline $fh) or return;
        defined(my $quality = scalar readline $fh) or return;
        chomp($readid_again, $quality);

        $readid =~ s/\^@//;
        $readid_again =~ s/\^\+//;
        return [$readid, $sequence, $readid_again, $quality];
    }
    elsif ($self->linesper()==2){
        $readid =~ s/\^>//;
        return [$readid, $sequence];
    }
    else{
        croak "only linesper 2,4 supported";
    }
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
