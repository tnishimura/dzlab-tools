package GFF::Parser::XS;
use Moose;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use GFF;
use Text::CSV_XS;
use autodie;
extends 'Parser';

has normalize => (
    is => 'ro',
    default => 1,
);

has csv => (
    is => 'rw',
    init_arg => undef,
);

sub BUILD{
    # warn "GFF::Parser::XS constructor called";
    my ($self) = @_;

    $self->csv(
        Text::CSV_XS->new({
                eol      => "\n",
                sep_char => "\t",
            }));
}

=head2 $p->next()

return the next gff record. returns undef when eof. If skip = 0 in constructor,
then return unparseable/comment lines as undef and pragmas as strings.

=cut 

sub next{
    my ($self) = @_;
    GFFLOOP:
    while (my $row = $self->csv->getline ($self->filehandle())) {
        if (scalar @$row != 9){
            next GFFLOOP;
        }
        elsif (defined $row->[0] && $row->[0] =~ /^\s*#/){
            next GFFLOOP;
        }
        else {
            return GFF->new(
                sequence         => $row->[0] eq q{} || $row->[0] eq q{.} ? undef : $row->[0],
                source           => $row->[1] eq q{} || $row->[1] eq q{.} ? undef : $row->[1],
                feature          => $row->[2] eq q{} || $row->[2] eq q{.} ? undef : $row->[2],
                start            => $row->[3] eq q{} || $row->[3] eq q{.} ? undef : $row->[3],
                end              => $row->[4] eq q{} || $row->[4] eq q{.} ? undef : $row->[4],
                score            => $row->[5] eq q{} || $row->[5] eq q{.} ? undef : $row->[5],
                strand           => $row->[6] eq q{} || $row->[6] eq q{.} ? undef : $row->[6],
                frame            => $row->[7] eq q{} || $row->[7] eq q{.} ? undef : $row->[7],
                attribute_string => $row->[8] eq q{} || $row->[8] eq q{.} ? undef : $row->[8],
            );
        }
    }
    return;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;


=head1 NAME
 
GFF::Parser - Parse a GFF file line by line
 
=head1 SYNOPSIS
 
    use GFF::Parser;
    my $p = GFF::Parser->new(file => "file.gff");
    while (my $gff = $p->next()){
        # ...
    }
  
=head1 DESCRIPTION
 
Parser for GFF.

=head1 SUBROUTINES/METHODS 

=over

=item GFF::Parser->new(file => 'file', skip => 1);

=item $gff = $p->next()


=back

=cut

