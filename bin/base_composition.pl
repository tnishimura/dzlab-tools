#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/../lib";
use GFF;
use GFF::Parser;
use FastaReader;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if ($opt_help || !($opt_gff xor $opt_window_size) || !$opt_input || !$opt_output);

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}

my $fr = FastaReader->new(file => $opt_input, slurp => 1);

if ($opt_gff){
    my $parser = GFF::Parser->new(skip => 1,file => $opt_gff);

    while (my $gff = $parser->next()){
        my $body = $fr->get($gff->sequence,$gff->start, $gff->end, coord => 'f', rc => $gff->strand eq '-' ? 1 : 0);
        my $len = $gff->end - $gff->start +1;
        say join "\t",
        $gff->sequence, $gff->source, $gff->feature, $gff->start, $gff->end, $gff->score // 0, $gff->strand // '+', $gff->frame // 0,
        count_chunk($body,$len) . (defined $gff->attribute_string ? ";" . $gff->attribute_string : "");
    }
}
else{
    for my $seqname (sort $fr->sequence_list()) {
        my $len = $fr->get_length($seqname);
        my $position = 1;
        INNER:
        while ($position <= $len){
            my $chunk = $fr->get($seqname, $position, $position+$opt_window_size -1, lenient => 1);
            my $end = $position + length($chunk) -1;
            say join "\t", $seqname, qw/. bc/, $position, $end, qw/. . ./, count_chunk($chunk, $opt_window_size);
            $position = $end+1;
        }
    }
}

if ($opt_output ne '-'){
    close \*STDOUT;
}

sub count_chunk{
    my ($chunk, $len) = @_;
    my $a = ($chunk=~tr/aA/aA/); my $a_ratio = $a / $len;
    my $c = ($chunk=~tr/cC/cC/); my $c_ratio = $c / $len;
    my $g = ($chunk=~tr/gG/gG/); my $g_ratio = $g / $len;
    my $t = ($chunk=~tr/tT/tT/); my $t_ratio = $t / $len;
    return sprintf("a_ratio=%.5f;c_ratio=%.5f;g_ratio=%.5f;t_ratio=%.5f;a=%d;c=%d;g=%d;t=%d",$a_ratio,$c_ratio,$g_ratio,$t_ratio,$a,$c,$g,$t); 
}

=head1 NAME

base_composition.pl - Given a gff annotation file or a window size, and a
corresponding reference fasta file, add the approrpiate ratios for each base to
the attribute field (column 9) of the gff and spit it back out.

=head1 SYNOPSIS

Count base composition by annotation:

 base_composition.pl -g annotation.gff -o annotation-with-count.gff reference.fasta 

Count base composition by fixed window size:

 base_composition.pl -w 50 -o annotation-with-count.gff reference.fasta 

=head1 OPTIONS

=over

=item  -g <file> | --gff <file>

GFF annotation file.  

=for Euclid
    file.type:        readable

=item  -w <size> | --window-size <size>

=item  <input> 

Input fasta file.

=for Euclid
    input.type:        readable

=item  -o <file> | --output <file>

output. Use '-' (without quotes) to print to screen.

=item --help | -h

=back

=cut




