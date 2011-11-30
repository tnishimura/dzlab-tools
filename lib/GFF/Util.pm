package GFF::Util;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;
use GFF;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(format_gff parse_gff is_gff);

=head2 parse_gff($line)

read a line, and return a gff hashref.  maps '.' dot columns and non-existenct attributes to undef. returns 0 on
non-parseable lines or comments. *Returns pragmas as strings*.  Check that ref $result eq 'HASH'

it's a member method so it can be inherited

=cut

sub parse_gff{
    my ($line, $normalize) = @_;

    $normalize //= 1;
     
    return 0 unless $line;
    $line =~ tr/\n\r//d;

    if ($line =~ m/^\s*#(#)?/o){
        return defined $1 ? $' : 0; # $' = regex postmatch
    }

    my @split = split /\t/, $line, 9;
    (carp "unparseable GFF line: $line" && return 0) unless @split == 9;

    my %accum;
    @accum{qw/sequence source feature start end score strand frame attribute_string/}
    = map { ($_ eq q{.} || $_ eq q{} ) ? undef : $_ } @split;

    if (defined $accum{sequence} && $normalize){
        if ($normalize == -1){
            $accum{sequence} = lc $accum{sequence};
        }
        else{
            $accum{sequence} = uc $accum{sequence};
        }
    }

    return GFF->new(%accum);
}

sub is_gff{
    my ($gff) = @_;
    return ref $gff eq 'GFF';
}

=head2 format_gff

 format_gff seq => 'sequence_name', start => 123, end => 456;

=cut

sub format_gff{
    croak "format_gff requires even number of args" unless @_%2 == 0;
    my %parts = @_;
    if (defined $parts{seq} && defined $parts{sequence}){
        croak "format_gff: can't specify both seq and sequence";
    }
    if (defined $parts{attr} && defined $parts{attribute}){
        croak "format_gff: can't specify both attr and attribute";
    }
    my $sequence  = delete $parts{sequence} // delete $parts{seq} // '.';
    my $source    = delete $parts{source} // '.';
    my $feature   = delete $parts{feature} // '.';
    my $start     = delete $parts{start} // '.';
    my $end       = delete $parts{end} // '.';
    my $score     = delete $parts{score} // '.';
    my $strand    = delete $parts{strand} // '.';
    my $frame     = delete $parts{frame} // '.';
    my $attribute = delete $parts{attribute} // delete $parts{attr} // '.';
    if (%parts){
        croak "unknown gff columns " . join ", ", sort keys %parts;
    }

    return join "\t",
    $sequence, 
    $source,   
    $feature,  
    $start,    
    $end,      
    $score,    
    $strand,   
    $frame,    
    $attribute,
}

1;
