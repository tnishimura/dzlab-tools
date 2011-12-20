#!/usr/bin/env perl
package GFF::Parser;
use Moose;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use GFF;
use GFF::Util;
use autodie;
use Scalar::Util qw/looks_like_number/;

has filename_or_handle => (
    is => 'ro',
    required => 1,
    init_arg => 'file',
);

has skip => (
    is => 'ro',
    default => 1
);

has done => (
    is => 'rw',
    default => 0
);

has normalize => (
    is => 'ro',
    default => 1,
);

has launder => (
    is => 'ro',
    default => 0,
);

has putback => (
    traits  => ['Array'],
    is      => 'ro',
    isa     => 'ArrayRef[GFF]',
    default => sub { [] },
    handles => {
        push_putback => 'push',
        shift_putback  => 'shift',
        has_putback  => 'count',
    },
);

# privates

has filehandle => (
    is => 'rw',
    init_arg => undef,
);

sub BUILD{
    my ($self) = @_;
    if (ref $self->filename_or_handle eq 'GLOB'){
        $self->filehandle($self->filename_or_handle);
    }
    elsif (!ref $self->filename_or_handle && -f $self->filename_or_handle ){
        open my $fh, '<', $self->filename_or_handle
            or croak "cannot open $self->filename_or_handle";
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

=head2 $p->next()

return the next gff record. returns undef when eof. If skip = 0 in constructor,
then return unparseable/comment lines as undef and pragmas as strings.

=cut 

sub next{
    my ($self) = @_;
    if ($self->has_putback){
        return $self->shift_putback;
    }
    my $launder = $self->launder();
    while (defined (my $line = scalar readline $self->filehandle)){
        my $gff = parse_gff($line, $self->normalize);
        if (!$self->skip || is_gff($gff)){
            if (
                (
                    $launder &&
                    do {
                        my $start = $gff->start; 
                        my $end = $gff->end;
                        looks_like_number($start) && looks_like_number($end) && $start <= $end
                    }
                )
                || 
                ! $launder
            ){
                return $gff;
            }
        }
    }
    $self->done(1);
    return;
}

sub rewind{
    my ($self,$gff) = @_;
    $self->push_putback($gff);
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

