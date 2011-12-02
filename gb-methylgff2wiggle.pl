#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use Launch qw/drain/;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use File::Temp qw/tempfile/;
use File::Path qw/make_path/;
use File::Spec::Functions qw/rel2abs catdir catfile/;
use File::Basename qw/basename dirname/;
use DZUtil qw/open_cached close_cached_all clean_basename/;
use GFF::Statistics qw/gff_detect_width/;

END {close STDOUT}
$| = 1;

sub capped_log10 { 
    my $n = shift; 
    my $x = $n == 0 ? 0 : log($n)/log(10); 
    return $x >= 5 ? 5 : $x;
}

sub methyl2wiggle{
    my ($file, $dir) = @_;
    my %wigs;
    my %done;
    my $detected_width = gff_detect_width $file;

    my $p = GFF::Parser->new(file => $file, normalize => -1);

    while (defined(my $gff = $p->next())){
        my ($seq, $start, $end, $score, $c, $t) 
        = ($gff->sequence(), $gff->start(), $gff->end(), $gff->score(), $gff->get_column('c') // 0, $gff->get_column('t') // 0,);

        my $width = $end - $start + 1;

        # okay for a single window to be different size b/c the last one may be smaller
        if ($width != $detected_width){
            if (exists $done{$seq}){ die "uneven widths"; }
            else{ $done{$seq} = 1; }
        }

        $score //= "0.0000";
        if ($opt_ct){ $score = ($c + $t) == 0 ? 0 : $c / ($c + $t); }

        my $clean_basename = clean_basename($file);

        # create wigs if ! exist
        for my $type (qw/methyl coverage/) {
            if (!exists $wigs{$seq}{$type}){
                my (undef, $tempfile) = tempfile(catfile($dir, "wig-$type-$clean_basename-XXXXXXX"), UNLINK => ! $opt_debug); 
                my $fh = open_cached('>', $tempfile);
                say $fh "variableStep  chrom=$seq  span=$detected_width";
                $wigs{$seq}{$type} = $tempfile;
            }
        }

        say {open_cached '>', $wigs{$seq}{methyl}}   "$start\t$score";
        say {open_cached '>', $wigs{$seq}{coverage}} sprintf("%d\t%4f", $start, capped_log10($c + $t));
    }

    close_cached_all();

    return %wigs;
}

sub process_wiggle2gff3_output{
    if (@_ % 2 != 0) { die "arg error" }
    my %opt = @_;

    my ($file, $original, $dir, $seq, $type, $prefix) = @opt{qw/wigfile original dir seq type prefix/};
    $prefix //= clean_basename($original);

    my $contents = drain('wiggle2gff3', '--path', $dir, $file);

    # grab the GFF line
    my @gff_lines = grep { scalar(@$_) == 9 } map { [split /\t/, $_] } split /\n/, $contents;
    if (@gff_lines != 1){
        die "wiggle2gff3 output gff malformatted?";
    }
    my @fields = @{pop @gff_lines};

    $fields[2] = "$prefix-$seq-$type";
    $fields[8] =~ s/Name=000;/Name=$prefix;/;
    return join "\t", @fields;
}

#######################################################################

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || ! defined $opt_input || ! defined $opt_tmp_dir;

if (defined $opt_feature_name){
    $opt_feature_name = clean_basename($opt_feature_name);
}

# stdout
if (defined $opt_output){
    open my $fh, '>', $opt_output;
    select $fh;
}

# make directory
my $opt_tmp_dir = rel2abs($opt_tmp_dir);
if (! -d $opt_tmp_dir){
    make_path($opt_tmp_dir);
}

# poof
my %text_wigs = methyl2wiggle($opt_input, $opt_tmp_dir); 

# convert each.
for my $seq (sort keys %text_wigs){
    for my $type (qw/methyl coverage/) {
        say process_wiggle2gff3_output(
            wigfile  => $text_wigs{$seq}{$type},
            original => $opt_input,
            dir      => $opt_tmp_dir,
            seq      => $seq,
            type     => $type,
            prefix   => $opt_feature_name,
        );
    }
}

=head1 SYNOPSIS

Usage examples:

 gb-methylgff2wiggle.pl -d output_dir -f featurename input.gff > output.wig

=head1 OPTIONS

=over

=item  <input>

=for Euclid
    input.type:        readable

=item  -o <filename> | --output <filename>

Default to STDOUT.

=item  --ct 

Calculate methylation score from 'c' and 't'  from column 9 instead of using
column 6.

=item  -d <dir> | --tmp-dir <dir>

Directory for intermediate and binary wig file.  

=item  -f <feature> | --feature-name <feature>

The feature name prefix in column 3 of output. Default to filename.

=item --debug

Keep intermediate wiggle text file.

=item --help | -h


=back

=cut
