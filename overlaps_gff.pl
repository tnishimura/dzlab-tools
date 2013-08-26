#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use lib "$FindBin::Bin/Tree-Range/lib";
use GFF::Tree;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use Scalar::Util qw/refaddr/;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || !$opt_gff || !$opt_input;

if ($opt_no_body && ! $opt_downstream && ! $opt_upstream){
    say STDERR "You specified -nb/--no-body, but not -down or -up. You have to look for something." ;
    exit 1;
}
if ($opt_proportion_threshold != 0 && ($opt_no_body || $opt_downstream || $opt_upstream)){
    say STDERR "Sorry, -p / --proportion-threshold is not currently compatible with --no-body, -down or -up." ;
    exit 1;
}

if ($opt_invert && ($opt_proportion_threshold != 0 || $opt_no_body || $opt_downstream || $opt_upstream)){
    say STDERR "Sorry, --invert is not currently compatible with -p, --no-body, -down or -up." ;
    exit 1;
}

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}

my $gt = GFF::Tree->new(file => $opt_gff, lenient => 1);
my $p  = GFF::Parser->new(file => $opt_input, normalize => 0);

my %annotation_has_seq;

while (defined(my $gff = $p->next())){
    my ($seq, $start, $end) = ($gff->sequence, $gff->start, $gff->end);

    if (! exists $annotation_has_seq{uc $seq}){ 
        # first time seeing $seq
        $annotation_has_seq{uc $seq} = $gt->has_tree(uc $seq);
        if (! $annotation_has_seq{uc $seq}){
            say STDERR "warning: annotation does not have any entries for $seq";
        }
    }

    my $outputted;

    if ($annotation_has_seq{uc $seq}){
        # shove them in a hash first, with address as keys so values are unique
        my %results_set;
        if (! $opt_no_body){
            $results_set{refaddr $_} = $_ for $gt->search_overlap($seq,$start,$end);
        }
        if ($opt_downstream > 0){
            $results_set{refaddr $_} = $_ for $gt->search_overlap($seq,$end + 1,$end + $opt_downstream);
        }
        if ($opt_upstream > 0){
            $results_set{refaddr $_} = $_ for $gt->search_overlap($seq,$start -  $opt_upstream, $start - 1);
        }

        if ($opt_invert){
            if (keys %results_set == 0){
                say "$gff";
            }
        }
        # commit note: following block not changed for adding invert block below, only indented
        else{
            my @results = values %results_set;

            for my $result (@results){
                my $overlap = $result->{overlap};
                if ($opt_proportion_threshold == 0 || 
                    ($opt_proportion_threshold == 1 && $overlap == $gff->length) || 
                    $overlap / $gff->length() >= $opt_proportion_threshold
                ){

                    if (defined $opt_tag && defined(my $locus = $result->{item}->get_column($opt_tag))){
                        if ($opt_tag_only){
                            $gff->attribute_string($locus);
                            say $gff;
                        }
                        else{
                            say "$gff;$opt_tag=$locus";
                        }
                        $outputted = 1;
                    }
                    else{
                        say STDERR "Warning: overlapping but no $opt_tag for " . $result->{item}; 
                    }
                }
                if ($outputted && $opt_one_each){
                    last;
                }
            }
        }

    }

    if (! $outputted && $opt_no_skip){
        #$gff->score(undef); # don't remember why this was desired...
        $gff->attribute_string(undef) if $opt_tag_only;
        say $gff;
    }
}


=head1 NAME

overlaps_gff.pl - for each input gff line, if it overlaps with an entry in the annotation, 
print the gff line with additional annotation information.  

=head1 SYNOPSIS

Usage examples:

 overlaps_gff.pl -g annotation.gff -k -t 'ID' -o output.gff input.gff

=head1 OPTIONS

=over

=item  -g <file> | --gff <file>

Annotation file.

=for Euclid
    file.type:        readable

=item  <input>

=for Euclid
    input.type:        readable

=item  -k | --no-skip

Print line from input.gff even if no overlap is found.  

=item  -t <tag> | --tag <tag>

Locus tag in --gff annotation file. Optional. Default to giving you the entire attribute string (col 9) from the annotation.

=item -o <file> | --output <file>

=for Euclid
    file.default:     '-'

=item  -p <thresh> | --proportion-threshold <thresh>

Input GFF lines needs at least this proportion covered by annotation GFF. For example,
if this is set to .50, and the input GFF is windowed by 50, at least 25 bases need to be overlapping
with an annotation in order to be reported.  If not, skipped (unless --no-skip is used, in which case it's 
reported with a null score). Default to 0, mean any overlap is enough.

=for Euclid
    thresh.default:     0
    thresh.type:        number, thresh >= 0 && thresh <= 1
    thresh.type.error:  <thresh> must be between 0 and 1

=item  -to  | --tag-only 

=item  --one-each

Only print one entry for multiply overlapping input.

=item  -up <bp> | --upstream <bp>

Search <bp> upstream of input windows.

=for Euclid
    bp.default:     0
    bp.type:        int, bp >= 0 
    bp.type.error:  <bp> must be greater than 0

=item  -down <bp> | --downstream <bp>

Search <bp> downstream of input windows.

=for Euclid
    bp.default:     0
    bp.type:        int, bp >= 0 
    bp.type.error:  <bp> must be greater than 0

=item --no-body | -nb

Do not look for overlaps with body only.  This is useful when looking for overlaps with 
upstream/downstream regions only.  For example, -nb -down 100 -up 100 searches for overlaps
with 100 bp upstream and downstream of a input window, but NOT the window itself.

=item --col8-marker | -m

Replace column 8 (normally the GFF 'frame' column) with 1 when input window has
an overlap.  This is useful when used with --no-skip, so that you can quickly
differentiation between windows with or without overlaps.

=item --invert | -v 

Filter out all gff entries which DO overlap with something in the annotation.

=item --help | -h

=back

=cut



