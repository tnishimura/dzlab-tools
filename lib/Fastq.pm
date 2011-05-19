package Fastq;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use Carp;
use autodie;    
use DBI;

has file => ( required => 1, is => 'ro');

# privates
has index => ( init_arg => undef, is => 'rw' );
has dbh => ( init_arg => undef, is => 'rw' );
has sth => ( init_arg => undef, is => 'rw' );

sub BUILD{
    my ($self) = @_;
    $self->index($self->file . ".index");

    if (! -f $self->file || ! -f $self->index){
        croak "$self->file or $self->index not readable?";
    }
    my $i = $self->index;

    my $dbh = DBI->connect("dbi:SQLite:dbname=$i","","", {RaiseError => 1, AutoCommit => 1});
    $dbh->do("PRAGMA cache_size = 80000");

    my $sth = $dbh->prepare("select sequence from fastq where id = ?");

    $self->dbh($dbh);
    $self->sth($sth);
}

sub get{
    my ($self, $id) = @_;
    $self->sth->execute($id);

    my ($sequence) = $self->sth->fetchrow_array;

    if (! defined $sequence){
        croak "couldn't fetch $id";
    }

    return $sequence;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;

