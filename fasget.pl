#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use FastaReader;

use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if ($opt_help || ($opt_coord ne 'f' && $opt_coord ne 'r'));

$opt_reverse_compliment //= 0;

if ($opt_output ne '-' ){
    open my $fh, '>', $opt_output;
    select $fh; 
}

my $r = FastaReader->new(file => $opt_fasta, slurp => ! $opt_no_slurp);

if (!defined $opt_sequence && ! defined $opt_range{start} && ! defined $opt_range{end}){
    my %lengths = $r->sequence_lengths;
    
    my $total = 0;
    for (sort keys %lengths){
        if ($opt_gff){
            say join "\t", $_, '.', 'chromosome', 1, $lengths{$_}, qw/. . ./, "Name=$_";
        }
        else{
            say "$_:\t" . $lengths{$_};
        }
        $total += $lengths{$_};
    }
    if (! $opt_gff){
        say "total:\t" . $total; 
    }
}
elsif (defined $opt_sequence && ! defined $opt_range{start} && ! defined $opt_range{end}){
    say $r->get_pretty($opt_sequence, $opt_sequence, undef, undef, rc    => $opt_reverse_compliment);
}
else{
    say $r->get_pretty(
        "$opt_sequence\_$opt_range{start}\_$opt_range{end}", 
        $opt_sequence, $opt_range{start}, $opt_range{end}, 
        coord => $opt_coord,
        rc    => $opt_reverse_compliment,
        base  => $opt_base);
}


=head1 NAME
 
fasget.pl
 
=head1 SYNOPSIS

 fasget.pl

=head1 REQUIRED ARGUMENTS

=over

=item  <fasta> 

=for Euclid
    fasta.type:        readable

=back

=head1 OPTIONS

=over

=item  -s <seq> | --sequence <seq>

=item  -r <start> <end> | --range <start> <end>

=for Euclid
    start.type:        int, start >= 1 
    start.type.error:  <fasta> must be > 1
    end.type:        int, end >= 1 
    end.type.error:  <fasta> must be > 1

=item  -c <coordinate> | --coord <coordinate>

Coordinate. Can be 'f' or 'r'. Default 'f'.

=for Euclid
    coordinate.default:     'f'

=item -[no-]rc | --reverse-compliment

Whether to reverse-complement chunk. Default off.

=for Euclid
    false: -no-rc

=item  -b <base> | --base <base>

The coordinate of the first base. Default 1.

=for Euclid
    base.default:     1
    base.type:        int, base >= 0 && base <= 1

=item -l | --no-slurp

Do not slurp entire file into memory. 


=item -o <file> | --output <file>

=for Euclid
    file.default:     '-'

=item  -g | --gff 

If used without any other arguments, prints out sequence sizes in GFF format.

=item -h | --help

=back

=cut

