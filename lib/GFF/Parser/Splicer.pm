# faster alternative interface to gff parsing

package GFF::Parser::Splicer;
use Moose;
use strict;
use warnings;
use Data::Dumper;
use Carp;
use autodie;
use 5.010_000;
use v5.14;
use IO::Handle;

extends 'Parser';

has normalize => (
    is => 'ro',
    default => 1,
);

has columns => (
    traits  => ['Array'],
    is      => 'ro',
    isa     => 'ArrayRef[Str]',
    default => sub { [] },
    handles => {
        all_columns     => 'elements',
    },
);

has colspecs => (
    traits  => ['Array'],
    is      => 'ro',
    default => sub { [] },
    handles => {
        all_colspecs     => 'elements',
        add_colspec     => 'push',
        count_colspecs   => 'count',
    },
);

sub BUILD{
    # warn "GFF::Parser::XS constructor called";
    my ($self) = @_;

    for (map { lc } $self->all_columns()) {
        $_ = lc $_;
        when ([qw/sequence seq/]) { $self->add_colspec(0) }
        when ([qw/source/])       { $self->add_colspec(1) }
        when ([qw/feature/])      { $self->add_colspec(2) }
        when ([qw/start/])        { $self->add_colspec(3) }
        when ([qw/end/])          { $self->add_colspec(4) }
        when ([qw/score/])        { $self->add_colspec(5) }
        when ([qw/strand/])       { $self->add_colspec(6) }
        when ([qw/frame/])        { $self->add_colspec(7) }
        when ([qw/attr attribute_string col9/]) { 
            $self->add_colspec(8) 
        }
        default {
            $self->add_colspec(qr/\Q$_\E=([^;]+)/);
        }
    }
}

=head2 $p->next()

return the next gff record. returns undef when eof. If skip = 0 in constructor,
then return unparseable/comment lines as undef and pragmas as strings.

=cut 

sub next{
    my ($self) = @_;
    my @colspecs = $self->all_colspecs;
    GFFLOOP:
    while (defined(my $line = $self->filehandle()->getline())) {
        chomp $line;
        if ($line =~ /^\s*#/){
            next GFFLOOP;
        }

        my @splice = (undef) x $self->count_colspecs;
        my @split = split /\t/, $line, 9;

        if (scalar @split != 9){
            next GFFLOOP;
        }

        my $index = 0;
        for my $spec (@colspecs) {
            if (ref $spec eq q{}){
                my $v = $split[$spec];
                if ($v ne '.' && $v ne ''){
                    if ($spec == 0 and $self->normalize()){
                        $v = uc $v;
                    }
                    $splice[$index] = $v;
                }
            }
            else { # $spec is qr//
                if ($split[8] =~ $spec){
                    $splice[$index] = $1;
                }
            }
            ++$index;
        }
        return \@splice;
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

