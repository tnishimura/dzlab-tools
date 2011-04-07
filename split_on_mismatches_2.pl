#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use autodie;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use feature 'say';
use File::Basename;
use File::Spec::Functions;

pod2usage(-verbose => 99,-sections => [('NAME', 'SYNOPSIS', 'OPTIONS', 'REQUIRED ARGUMENTS')] )
    if $opt_help || ! ($opt_label xor ($opt_output_a && $opt_output_b));

open my $ain, '<', $opt_input_a;
open my $bin, '<', $opt_input_b;

my ($basename_a, $path_a, undef) = fileparse($opt_input_a,".sorted");
my ($basename_b, $path_b, undef) = fileparse($opt_input_b,".sorted");

my $output_file_a = $opt_output_a // catfile($path_a, $basename_a) . ".$opt_label.filtered";
my $output_file_b = $opt_output_b // catfile($path_b, $basename_b) . ".$opt_label.filtered";

open my $aout, '>', $output_file_a;
open my $bout, '>', $output_file_b;

open my $aerror, '>', "$output_file_a.$opt_error_suffix";
open my $berror, '>', "$output_file_b.$opt_error_suffix";

CMP:
while (defined (my $a_record = <$ain>) and 
    defined (my $b_record = <$bin>)) {
    $a_record =~ tr/\n\r//d;
    $b_record =~ tr/\n\r//d;

    my ($a_id, $a_strand, $a_mm, $a_rawcoord) = (split /\t/, $a_record)[0..3];
    my ($b_id, $b_strand, $b_mm, $b_rawcoord) = (split /\t/, $b_record)[0..3];

    die "$opt_input_a or $opt_input_b not sorted??" if ($a_id ne $b_id); 

    $a_mm = get_score ($a_mm);
    $b_mm = get_score ($b_mm);

    #say $a_rawcoord;
    #say $b_rawcoord;

    if (! defined $a_mm and ! defined $b_mm) {
        next CMP; # read matched nothing in either ecotypes
    }
    elsif (defined $a_mm && defined $b_mm ){
        if ($a_rawcoord =~ s/.*chr\w:(\d+).*/$1/xmsi &&
            $b_rawcoord =~ s/.*chr\w:(\d+).*/$1/ixms){

            # both mapped somewhere
            if ($a_rawcoord == $b_rawcoord && $a_strand eq $b_strand){
                # both mapped to same place
                next if $a_mm == $b_mm; # tied.
                say $aout $a_record if $a_mm < $b_mm; # A wins
                say $bout $b_record if $a_mm > $b_mm; # B wins
            } else {
                # mapped to different places
                say $aerror $a_record;
                say $berror $b_record;
            }
        } else{
            # hackish-- reaching this means that rawcoord must've been negative?
            # fix in parse_bowtie
            say $aerror $a_record;
            say $berror $b_record;
        }
    }
    elsif (! defined $a_mm) {
        # did match somewhere for b but not a
        say $bout $b_record;
    }	
    elsif (! defined $b_mm) {
        # did match somewhere for a but not b
        say $aout $a_record;
    }
    else {croak "Impossible situation:\n$a_record\n$b_record"}
}

close $aout; close $bout;
close $ain;   close $bin;
close $aerror;
close $berror;

sub get_score {
    my ($mm) = @_;

    return if 'NM' eq $mm;
    my @mm = split /:/, $mm;

    for my $i (0 .. @mm - 1) {
        return $i if 1 == $mm[$i];
    }
}

__END__


=head1 NAME

split_on_mismatches_2.pl - Filter sorted bowtie inputs into two files based on mismatch counts.

=head1 SYNOPSIS

If outputs not explicitly given, Produces four files:
 
 left.eland.<label>.filtered
 left.eland.<label>.filtered.error
 right.eland.<label>.filtered
 right.eland.<label>.filtered.error

overwrites input files with correct imprinting alignments:

 perl split_on_mismatches_2.pl -l CxL -a left.eland -b right.eland 

=head1 REQUIRED ARGUMENTS

=over

=item  -a <eland> | --input-a <eland>

Ecotype A eland file

=for Euclid eland.type:        readable

=item  -b <eland> | --input-b <eland>

Ecotype B eland file

=for Euclid eland.type:        readable


=back

=head1 OPTIONS

=over

=item  -e <suffix> | --error-suffix <suffix>

=for Euclid
    suffix.default:     'error'

=item  -oa <output> | --output-a <output>

Filtered A eland output file

=item  -ob <output> | --output-b <output>

Filtered B eland output file

=item  -l <label> | --label <label>

Label for split.

=item --help | -h

=back

=cut

