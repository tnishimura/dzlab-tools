package Eland::Parser;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Moose;
use Carp;
use autodie;    

with 'ParserRole';

has 'fastareader' => (
    is => 'ro',
);

sub next{
    my $self = shift;
    if (defined (my $line = scalar readline $self->filehandle)){
        return parse_eland($line, $self->fastareader());
    }
    return;
}

sub parse_eland{
    my ($line, $fasta_reader) = @_;
    chomp $line;

    my ($read_id, $read_sequence, my $match_counts_string, my $matches_string) = split /\t/, $line;

    my @positions;

    if ($match_counts_string ne 'NM'){
        for my $match (split /,/, $matches_string) {
            my $reverse = $match =~ s/^RC_//;
            if ($match =~ m/([^:]+):(\d+)(F|R)(\d)/){
                my $seqid = $1;

                my $coord_start = $2;
                my $coord_end = $coord_start + length($read_sequence) - 1;

                my $mismatch = $4;

                if (defined $fasta_reader){
                    if ($reverse){
                        ($coord_end, $coord_start) = ($fasta_reader->reverse2forward($seqid, $coord_start), $fasta_reader->reverse2forward($seqid, $coord_end));
                    }
                    push @positions, [$seqid, $mismatch, $reverse, $coord_start, $coord_end];
                }
                else{
                    push @positions, [$seqid, $mismatch];
                }
            }
            else{
                die "$line unparseable";
            }
        }
    }
    return [$read_id, $read_sequence, @positions];
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
