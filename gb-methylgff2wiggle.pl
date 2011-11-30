#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use Launch qw/cast/;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use File::Temp qw/tempfile tempdir/;
use File::Path qw/make_path/;
use File::Spec::Functions qw/rel2abs catdir catfile/;
use File::Basename qw/basename dirname/;

END {close STDOUT}
$| = 1;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) if $opt_help || ! defined $opt_input;

#######################################################################
# temp dir for compiled wiggle-- either given my --tmp-dir or $HOME/.widdgle

my $tempdir;
if (defined $opt_tmp_dir){
    if (! -d $opt_tmp_dir){
        make_path($opt_tmp_dir);
    }
    $tempdir = rel2abs($opt_tmp_dir)
}
else{
    $tempdir = catdir($ENV{HOME}, ".wiggle");
    if (! -d $tempdir){
        make_path($tempdir);
    }
}

#my $track  = $opt_trackname // 'track';
my $p      = GFF::Parser->new(file => $opt_input, normalize => -1);

#######################################################################
# process gff line by line

my %seen;
my %done;

while (defined(my $gff = $p->next())){
    my ($seq, $start, $end, $score, $c, $t) 
    = ($gff->sequence(), $gff->start(), $gff->end(), $gff->score(), $gff->get_column('c'), $gff->get_column('t'),);

    my $width = $end - $start + 1;

    # okay for a single window to be different size b/c the last one may be smaller
    if ($width != $opt_width){
        if (exists $done{$seq}){ die "uneven widths"; }
        else{ $done{$seq} = 1; }
    }

    $score //= "0.0000";
    if ($opt_ct){ $score = $c / ($c + $t); }
    
    # initialize wig
    if (!exists $seen{$seq}){
        my (undef, $tempfile) = 

        # if debug on, then in same dir as compiled wigs.  otherwise in /tmp
        $opt_debug ? tempfile(catfile($tempdir, "wig-$seq-XXXXXXX"), UNLINK => 0) : tempfile(UNLINK => 1);
        open my $fh, '>', $tempfile;

        say STDERR "creating $tempfile for $seq";

        say $fh "variableStep  chrom=$seq  span=$opt_width";
        $seen{$seq} = {fh => $fh, filename => $tempfile};
    }

    my $output_fh = $seen{$seq}{fh};
    say $output_fh "$start\t$score";
}

#######################################################################
# wiggle2gff3 for each sequence, plus post-processing.

# create a temporary output file
# if debug on, then in same dir as compiled wigs.  otherwise in /tmp
my (undef, $tmpout) = $opt_debug ? tempfile(catfile($tempdir, "gff-XXXXXXX"), UNLINK => 0) : tempfile(UNLINK => 1);

# convert each.
for my $seq (sort keys %seen){
    my ($fh, $filename) = @{$seen{$seq}}{qw/fh filename/};
    close $fh;

    cast("wiggle2gff3 --path $tempdir $filename >> $tmpout");
}

# eliminate all non-gff lines (like comments) in intermediate gff file:
cast("perl -i -wlnaF'\\t' -e '\@F == 9 and print' $tmpout");

# change column 3
if (! defined $opt_feature_name){
    $opt_feature_name = basename($opt_input, qw/.gff .GFF/);
}

cast("perl -i -wlpe 's/\tmicroarray_oligo\t/\t$opt_feature_name\t/' $tmpout");

#######################################################################
# output gff

# slurp
my $content = do {
    local $/;
    open my $tmpfh, '<', $tmpout;
    my $content = scalar <$tmpfh>;
    close $tmpfh;
    $content;
};

# spit
if (defined $opt_output){
    open my $outfh, '>', $opt_output;
    select $outfh;
}
print $content;

=head1 NAME

methylgff2wiggle.pl - ...

=head1 SYNOPSIS

Usage examples:

 methylgff2wiggle.pl [options]...

=head1 OPTIONS

=over

=item  -i <filename> | --input <filename>

=for Euclid
    filename.type:        readable

=item  -o <filename> | --output <filename>

=item  --ct 

Calculate methylation score from 'c' and 't'  from column 9 instead of using
column 6.

=item  -d <dir> | --tmp-dir <dir>

Directory for binary wig file.  Also, intermediate files if --debug on

=item  -f <feature> | --feature-name <feature>

The feature name in column 3 of output. Default to filename.

=item  -t <name> | --trackname <name>

=item  -w <size> | --width <size>

=for Euclid
    size.default:     50

=item  --debug 

Put intermediate files in --tmp-dir

=item --help | -h

=back

=cut
