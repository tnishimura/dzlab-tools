#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Pod::Usage;
use Getopt::Euclid qw( :vars<opt_> );
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;


pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
unless ($opt_annotation && scalar %opt_files);

my @files = @opt_files{qw/f1 f2 f3 f4/};

my @labels = scalar %opt_labels ? @opt_labels{qw/f1 f2 f3 f4/} : map {basename $_} @files;

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}

my %table;
# iinitailly fill with zero's
{
    my $p = GFF::Parser->new(file => $opt_annotation);
    while (defined(my $gff = $p->next())){
        my $locus = $gff->get_column($opt_locus_tag);
        if ($locus){
            $table{$locus} = [qw/0 0 0 0/];
        }
    }
}

my $col = 0;
for my $file (@files) {
    # body...
    my $p = GFF::Parser->new(file => $file);
    while (defined(my $gff = $p->next())){
        my $locus = $gff->get_column($opt_locus_tag);
        if ($locus){
            my $score = $gff->score // 0;
            $table{$locus}[$col] += $score;
        }
    }
    $col++;
}

say join "\t", "locus", @labels;
for my $locus (sort keys %table) {
    say join "\t", $locus, @{$table{$locus}};
}


=head1 NAME
 
divorce-gene-table.pl
 
=head1 SYNOPSIS

 divorce-gene-table.pl

=head1 OPTIONS

=over

=item -o <file> | --output <file>

=for Euclid
    file.default:     '-'

=item  -f <f1> <f2> <f3> <f4> | --files <f1> <f2> <f3> <f4> 

The files

=for Euclid
    f1.type:        readable
    f2.type:        readable
    f3.type:        readable
    f4.type:        readable

=item  -l <f1> <f2> <f3> <f4> | --labels <f1> <f2> <f3> <f4> 

Column labels. defaults

=item  -a <annotation> | --annotation <annotation>

GFF annotation

=for Euclid
    annotation.type:        readable

=item  -t <locus> | --locus-tag <locus>

Default 'ID'

=for Euclid
    locus.default:     'ID'

=item -h | --help

=back

=cut

