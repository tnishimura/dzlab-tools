package Eland::Parser;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(parse_eland);
our @EXPORT = qw();

sub parse_eland{
    my ($line, $fasta_reader) = @_;
    chomp $line;

    my %results;

    ($results{read_id}, $results{read_sequence}, my $match_counts_string, my $matches_string) = split /\t/, $line;
    $results{positions} = [ ];

    if ($match_counts_string eq 'NM'){
        $results{NM} = 1;
    }
    else{
        $results{NM} = 0;
        for my $match (split /,/, $matches_string) {
            my $reverse = $match =~ s/^RC_//;
            if ($match =~ m/([^:]+):(\d+)(F|R)(\d)/){
                my $seqid = $1;

                my $coord_start = $2;
                my $coord_end = $coord_start + length($results{read_sequence}) - 1;

                my $mismatch = $4;

                if ($reverse){
                    ($coord_end, $coord_start) = ($fasta_reader->reverse2forward($seqid, $coord_start), $fasta_reader->reverse2forward($seqid, $coord_end));
                }
                push @{$results{positions}}, [$seqid, $mismatch, $reverse, $coord_start, $coord_end];
            }
            else{
                die "$line unparseable";
            }
        }
    }
    return \%results;
}

1;

