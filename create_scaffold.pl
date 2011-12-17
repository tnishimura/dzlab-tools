#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF;
use GFF::Parser;
use FastaReader;

use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if ($opt_help || !$opt_gff || !$opt_input || !$opt_output);

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}

my $fr  = FastaReader->new(file => $opt_input, slurp => 1);
my $parser = GFF::Parser->new(skip => 1,file => $opt_gff);

my $strand_error = 0;
while (my $gff = $parser->next()){
    my $str = $gff->strand;
    if (defined $str && ($str ne '+' || $str ne '-') ){ # IGNORE entries without strands
        print $fr->get_pretty(
            ($opt_full ? 
                join "_", $gff->get_column($opt_locus_id), $gff->sequence, $gff->start, $gff->end 
                : $gff->get_column($opt_locus_id)
            ),
            $gff->sequence,
            $gff->start, 
            $gff->end, 
            coord => 'f', 
            rc => $str eq '-' ? 1 : 0
        );
    }
    else {
        $strand_error++;
    }
        #print format_fasta(
        #    #sprintf("%s | %s, %s, %s", $gff->get_column($opt_locus_id),$gff->start, $gff->end, $gff->strand),
        #    $gff->get_column($opt_locus_id),
        #    fasta_subseq($fasta,$gff->sequence,$gff->start, $gff->end, coord => 'f', rc => $gff->strand eq '-' ? 1 : 0)
        #);
}
if ($opt_output ne '-'){
    close \*STDOUT;
}
if ($strand_error){
    say STDERR "warning: ignored $strand_error lines with an undefined/malformed strand columns";
}

=head1 NAME

create_scaffold.pl - Create a scaffold given a fasta file and a gff annotation.

=head1 SYNOPSIS

Usage examples:

 create_scaffold.pl -g annotation.gff -o scaffold.fasta reference.fasta 

=head1 OPTIONS

=over

=item  -g <file> | --gff <file>

GFF annotation file.  

=for Euclid
    file.type:        readable

=item  -i <id> | --locus-id <id>

Name of the attribute in column 9 of GFF. Used to name the scaffolds. Default 'ID'.

=for Euclid
    id.default:     'ID'

=item  <input> 

Input fasta file.

=for Euclid
    input.type:        readable

=item  -o <file> | --output <file>

output. Use '-' (without quotes) to print to screen.

=item  -f | --full

=item --help | -h

=back

=cut




