#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config bundling);
use Pod::Usage;
use List::Util qw/max min sum/;

use Statistics::Descriptive;
#use Devel::Size qw(size total_size);

my $output;
my $width   = 50;
my $step    = 50;
my $merge   = 0;
my $no_sort = 0;
my $no_skip = 0;
my $scoring = 'meth';    # meth, average, sum or seq_freq
my $reverse = 0;         # reverse score and count
my $gff;
my $tag = 'ID';
my $qtag;
my $feature;
my $absolute = 0;

# Grabs and parses command line options
my $result = GetOptions(
    'width|w=i'    => \$width,
    'step|s=i'     => \$step,
    'scoring|c=s'  => \$scoring,
    'merge|m=s'    => \$merge,
    'no-sort|n'    => \$no_sort,
    'no-skip|k'    => \$no_skip,
    'gff|g=s'      => \$gff,
    'tag|t=s'      => \$tag,
    'collect-tag|u=s' => \$qtag,
    'feature|f=s'  => \$feature,
    'reverse|r'    => \$reverse,
    'absolute|b=s' => \$absolute,
    'output|o=s'   => \$output,
    'verbose'      => sub { use diagnostics; },
    'quiet'        => sub { no warnings; },
    'help'         => sub { pod2usage( -verbose => 1 ); },
    'manual'       => sub { pod2usage( -verbose => 2 ); }
);

my %scoring_dispatch = (
    meth                 => \&fractional_methylation,
    sum                  => \&sum_scores,
    average              => \&unweighed_average,
    weighed_average      => \&full_weighed_average,
    half_weighed_average => \&half_weighed_average,
    locus_collector      => \&locus_collector,
    microarray           => \&microarray_average,
);

# Check required command line parameters
pod2usage( -verbose => 99,-sections => [qw/NAME DESCRIPTION SYNOPSIS OPTIONS SCORING/]  )
    unless @ARGV
        and $result
        and ( ( $width and $step ) or $gff )
        and exists $scoring_dispatch{$scoring};

if ($output) {
    open my $USER_OUT, '>', $output
        or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

if ($merge) {
    $width = 1;
    $step  = 1;
}

if ($absolute) {
    croak
        '-b, --absolute option only works with arabidopsis or rice or puffer'
        unless (lc $absolute eq 'arabidopsis'
                or lc $absolute eq 'rice'
                or lc $absolute eq 'puffer');
    $no_skip = 1;
}

my $fasta_lengths = {};
if ( $absolute and $no_skip ) {
    $fasta_lengths = index_fasta_lengths( handle => 'DATA' );
}

my $window_iterator = sub { };
if ($gff) {
    $window_iterator = make_annotation_iterator( # really? annonation?
        file  => $gff,
        tag   => $tag,
        merge => $merge
    );
}

##########################################################
# get chromosome names from data file

open my $GFF, '<', $ARGV[0] or croak "Can't open $ARGV[0]:$!";
my $gff_iterator = make_gff_iterator( parser => \&gff_read, handle => $GFF );
my %chromosomes;

print STDERR "Loading groups...";
while ( my $gff_line = $gff_iterator->()) {
    next unless ref $gff_line eq 'HASH';    
    my $sequence = $gff_line->{seqname};

    ++$chromosomes{$sequence} and print STDERR "\n$sequence"
    unless exists $chromosomes{$sequence};
}
close $GFF or croak "Can't close $ARGV[0]:$!";
print STDERR "...done\n";


SEQUENCE:
for my $sequence ( sort keys %chromosomes ) {

    open my $GFF, '<', $ARGV[0] or croak "Can't open $ARGV[0]:$!";
    my $gff_iterator = make_gff_iterator( parser => \&gff_read, handle => $GFF );
    
    my %gff_records = ();

    ##########################################################
    # slurp GFF into % gff_records. should be factored out
    
    print STDERR "Loading $sequence...";
  LOAD:
    while ( my $gff_line = $gff_iterator->() ) {
	next LOAD unless ref $gff_line eq 'HASH';	
	if ($gff_line->{seqname} eq $sequence) {
	    push @{ $gff_records{ $sequence } },
	    { map{ $_ => $gff_line->{$_} } qw(start end score attribute) };  
            # does this need more? $seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute
	}
    }
    close $GFF or croak "Can't close $ARGV[0]:$!";
    print STDERR "done\n";

    #printf STDERR "gff_records occupies %g MiB of memory\n",
    #total_size( \%gff_records ) / 1024 / 1024;

    ##########################################################
    # sort each gff_record by start

    print STDERR 'Sorting ', scalar @{ $gff_records{$sequence} }, ' records...';
    unless ($no_sort) {
        @{ $gff_records{$sequence} }
            = sort { $a->{start} <=> $b->{start} }
            @{ $gff_records{$sequence} };
    }
    print STDERR "done\n";

    unless ($gff) {
        $window_iterator = make_window_iterator(
            width => $width,
            step  => $step,
            lower => 1,
            upper => (
                        $absolute
                    and $no_skip
                    and %$fasta_lengths
                    and exists $fasta_lengths->{$absolute}{$sequence}
                )
            ? $fasta_lengths->{$absolute}{$sequence}
            # this doesn't seem right... even if sorted, gff_records is sorted
            # by sequence start, not end. should we find the max end manually?
            : $gff_records{$sequence}[-1]->{end},
        );
    }

    print STDERR "Windowing...";
  WINDOW:
    while ( my ( $targets, $locus ) = $window_iterator->($sequence) ) {

        my $brs_iterator = binary_range_search(
            targets => $targets,
            queries => $gff_records{$sequence},
        );

        my $scores_ref = $scoring_dispatch{$scoring}->($brs_iterator, 
            reverse => $reverse, 
            targets => $targets,
            locus_tag => $qtag,
        );

        my ( $score, $attribute );

        if ( ref $scores_ref eq 'HASH' ) {

            $score = sprintf( "%g", $scores_ref->{score} );
            delete $scores_ref->{score};

            $attribute = "ID=$locus; " if $locus;
            $attribute .= join q{; },
                map { "$_=" . $scores_ref->{$_} }
                sort keys %{$scores_ref};

        }
        elsif ( ref $scores_ref eq 'ARRAY' ) {
            $score = q{.};
            $attribute = "ID=$locus; " if $locus;
            $attribute .= join q{; }, @$scores_ref;
        }
        elsif ($no_skip) {
            $score = q{.};
            $attribute = $locus ? "ID=$locus" : q{.};
        }
        else {
            next WINDOW;
        }

        print join( "\t",
            $sequence, 'dzlab',
            ( $feature ? $feature : $gff ? 'locus' : 'window' ),
            $targets->[0][0], $targets->[-1][1],
            $score, q{.}, q{.}, $attribute, ),
            "\n";
    }
    print STDERR "done\n";

    delete $gff_records{$sequence};    
}

##########################################################
# index_fasta_lengths (from __DATA__ below)
#
# return {species => {sequence => length}}

sub index_fasta_lengths {
    my %options = @_;

    my $handle = $options{handle};

    if ( $options{file} ) {
        open $handle, '<', $options{file} or croak $!;
    }

    my %fasta_lengths;
    my $active;
    while (<$handle>) {
        chomp;
        if (m/^\s*#/) {
            s/#//;
            $active = lc $_;
        }
        else {
            my ( $sequence, $length ) = split /\t/;
            $fasta_lengths{$active}{$sequence} = $length;
        }
    }

    if ( $options{file} ) {
        close $handle or croak $!;
    }

    return \%fasta_lengths;
}

##########################################################
# make_window_iterator
# iterate from lower to upper by step.

sub make_window_iterator {
    my (%options) = @_;

    my $width = $options{width} || croak 'Need window parameter';
    my $step  = $options{step}  || croak 'Need step parameter';
    my $lower = $options{lower} || croak 'Need lower bound parameter';
    my $upper = $options{upper} || croak 'Need upper bound parameter';

    return sub {
        my $i = $lower;
        $lower += $step;

        if ( $i <= $upper - $width + 1 ) {
            return [ [ $i, $i + $width - 1 ] ];
        }
        elsif ( $i < $upper ) {
            return [ [ $i, $upper ] ];
        }
        else {
            return;
        }
    }
}

##########################################################
# make_annotation_iterator
#
# return an iterator which, given a sequence, will return 
# ( [ ranges ], $locus ), on every call, until loci runs out
sub make_annotation_iterator {
    my (%options) = @_;

    my $annotation_file = $options{file};
    my $locus_tag       = $options{tag} || 'ID';
    my $merge_exons     = $options{merge} || 0;

    open my $GFFH, '<', $annotation_file
        or croak "Can't read file: $annotation_file";

    my $gff_iterator
        = make_gff_iterator( parser => \&gff_read, handle => $GFFH );

    my $annotation = index_annotation( iterator => $gff_iterator, %options );

    close $GFFH or croak "Can't close $annotation_file: $!";

    my %annotation_keys = ();
    return sub {
        my ($sequence) = @_;
        return unless $sequence;

        unless ( exists $annotation_keys{$sequence} ) {
            @{ $annotation_keys{$sequence} }
                = sort keys %{ $annotation->{$sequence} };
        }

        my $locus = shift @{ $annotation_keys{$sequence} };
        return unless $locus;
        my $ranges = $annotation->{$sequence}{$locus};
        delete $annotation->{$sequence}{$locus};

        if ( $locus and @$ranges ) {
            return [ uniq_ranges($ranges) ], $locus;
        }
        else {
            return;
        }
    };
}

##########################################################
# index annotations
# 
# accept a gff iterator, a locus tag ('ID') and a feature (column 3) to merge
# on.
#
# return { sequence_names => { locus_id => [ranges]} }

sub get_locus_id{
    my ($locus_gff, $locus_tag,$merge_feature, $first) = @_;

    my ($locus_id) = $locus_gff->{attribute} =~ m/$locus_tag[=\s]?([^;,]+)/;


    if ( !defined $locus_id ) {
        ( $locus_id, undef ) = split /;/, $locus_gff->{attribute} if $first;
        $locus_id ||= q{.};
    }
    else {
        $locus_id =~ s/["\t\r\n]//g;
        $locus_id
        =~ s/\.\w+$// # this fetches the parent ID in GFF gene models (eg. exon is Parent=ATG101010.n)
        if $merge_feature and $locus_gff->{feature} eq $merge_feature;
    }
    return $locus_id;
}

sub index_annotation {
    my (%options) = @_;

    my $gff_iterator = $options{iterator}
        || croak 'Need an iterator parameter';
    my $locus_tag     = $options{tag}   || 'ID';
    my $merge_feature = $options{merge} || 0;

    my %annotation = ();

LOCUS:
    while ( my $locus = $gff_iterator->() ) {

        next LOCUS unless ref $locus eq 'HASH';

        my $locus_id = get_locus_id($locus, $locus_tag, $merge_feature,1);

        push @{ $annotation{ $locus->{seqname} }{$locus_id} },
            [ $locus->{start}, $locus->{end} ]
            unless ( $merge_feature and $locus->{feature} ne $merge_feature );

    }

    return \%annotation;
}

##########################################################
# remove duplicate ranges
# 
# $ranges is an array ref of ranges
# [[$start_a, $end_a], [$start_b, $end_b], ... ]
sub uniq_ranges {
    my ($ranges) = @_;

    my %seen;
    my @uniq;

    for (@$ranges) {
        push @uniq, $_ unless $seen{"$_->[0];$_->[1]"}++;
    }

    return wantarray ? @uniq : [@uniq];
}

##########################################################
# for the ranges from brs, count the c's and the t's in the attributes field

sub fractional_methylation {
    my ($brs_iterator,%opt) = @_;

    my ( $c_count, $t_count, $score_count ) = ( 0, 0, 0 );
    my $lc = make_locus_collector($opt{locus_tag}) if $opt{locus_tag};

  COORD:
    while ( my $gff_line = $brs_iterator->() ) {

        next COORD unless ref $gff_line eq 'HASH';
        $lc->($gff_line) if $opt{locus_tag};

        my ( $c, $t ) = $gff_line->{attribute} =~ m/
                                                     c=(\d+)
                                                     .*
                                                     t=(\d+)
                                                 /xms;

        $c_count += $c;
        $t_count += $t;
        $score_count++;
    }

    if ( $score_count and ( $c_count or $t_count ) ) {
        return {
            score => $c_count / ( $c_count + $t_count ),
            c     => $c_count,
            t     => $t_count,
            n     => $score_count,
            ($opt{locus_tag} ? $lc->() : ())
        };
    }
}

##########################################################
# return sum of scores of ranges returned by brs

sub sum_scores {
    my ($brs_iterator, %opt) = @_;
    my $reverse = $opt{reverse};

    my ( $score_sum, $score_count ) = ( 0, 0 );

    my $lc = make_locus_collector($opt{locus_tag}) if $opt{locus_tag};

COORD:
    while ( my $gff_line = $brs_iterator->() ) {
        next COORD unless ref $gff_line eq 'HASH';
        $lc->($gff_line) if $opt{locus_tag};

        $score_sum += $gff_line->{score} eq q{.} ? 0 : $gff_line->{score};
        $score_count++;
    }

    if ($reverse) {
        my $tmp      = $score_sum;
        $score_sum   = $score_count;
        $score_count = $tmp;
    }

    if ($score_count or ($reverse and $score_sum)) {
        return {
            score => $score_sum,
            n     => $score_count,
            ($opt{locus_tag} ? $lc->() : ())
        };
    }
}

##########################################################
# seq_freq 
# from the ranges from brs_iterator, calculate the
# frequency of each base at each position. (?)

sub seq_freq {
    my ($brs_iterator) = @_;

    my ( $seq_freq, $seq_count, $score_avg ) = ( undef, 0, 0 );

  COORD:
    while ( my $gff_line = $brs_iterator->() ) {

        next COORD unless ref $gff_line eq 'HASH';

        $score_avg = $score_avg
            + ( $gff_line->{score} - $score_avg ) / ++$seq_count;

        my ($sequence) = $gff_line->{attribute} =~ m/seq=(\w+)/;
        my @sequence   = split //, $sequence;

        # for each position, tally the base
        for my $i (0 .. @sequence - 1) {

            for (qw/a c g t/) {
                $seq_freq->[$i]{$_} = 0 unless exists $seq_freq->[$i]{$_}
            }
        
            $seq_freq->[$i]{lc $sequence[$i]}++
        }
        $seq_count++;
    }

    # attribute = { position_id => {a => a_ratio, c => c_ratio, ...}}
    my $i = 1;
    my %attribute;
    for my $position (@$seq_freq) {

        my $total = sum values %$position;

        $attribute{$i++}
        = join q{,}, map {
            sprintf "%g", ($position->{$_} / $total)
        } sort keys %$position;
    }

    if ($seq_count and $seq_freq) {
        return {
            %attribute,
            score => $score_avg,
            n     => $seq_count,
        };
    }
}

#==================================================================
# Locus

sub locus_collector{
    my ($brs_iterator,%opt) = @_;
    my $locus_tag = $opt{locus_tag};
    my $counter = 0;

    my $lc = make_locus_collector($locus_tag);
    COORD:
    while ( my $gff_line = $brs_iterator->() ) {
        next COORD unless ref $gff_line eq 'HASH';
        ++$counter;
        $lc->($gff_line);
    }

    if ($counter) {
        return {
            score       => $counter,
            $lc->()
        };
    }
}
sub make_locus_collector{
    my $locus_tag = shift;
    my @accum;
    my $unknown = 0;
    return sub {
        if (@_){
            my $gff_line = shift;
            my $l = get_locus_id($gff_line, $locus_tag,0,0);
            if ('.' eq $l){ 
                ++$unknown;
            } else {
                push @accum, $l;
            }
        } else {
            return (
                ID => join(",", @accum),
                unknowns    => $unknown
            );
        }
    };
}

#==================================================================
# weight_average
# return size of overlap
sub overlap{
    my ($x_start, $x_end, $y_start, $y_end) = @_;

    # no overlap
    if ($y_end < $x_start || $x_end < $y_start) {return 0; } 
    # x: |------|
    # y:  |--|
    elsif ($x_start <= $y_start && $y_end <= $x_end){ return $y_end - $y_start + 1 } #complete overlap
    # x:  |--|
    # y: |------|
    elsif ($y_start <= $x_start && $x_end <= $y_end){ return $x_end - $x_start + 1 } #complete overlap
    # x:    |------|
    # y:  |---|
    elsif ($y_start <= $x_start && $y_end <= $x_end) { return $y_end-$x_start +1} #partial
    # x: |------|
    # y:     |---|
    elsif ($x_start <= $y_start && $x_end <= $y_end) { return $x_end-$y_start +1} #partial
    # shouldn't happen
    else { return 0 }
}

sub full_weighed_average{
    my ($brs_iterator, %opt) = @_;
    weighed_average($brs_iterator, %opt, source => 1, dest => 1);
}

sub half_weighed_average{
    my ($brs_iterator, %opt) = @_;
    weighed_average($brs_iterator, %opt, source => 1, dest => 0);
}
sub unweighed_average{
    my ($brs_iterator, %opt) = @_;
    weighed_average($brs_iterator, %opt, source => 0, dest => 0);
}

sub weighed_average{
    my ($brs_iterator, %opt) = @_;
    my $targets = $opt{targets};
    my $stat = Statistics::Descriptive::Full->new();
    my $counter = 0;
    my $lc = make_locus_collector($opt{locus_tag}) if $opt{locus_tag};

    COORD:
    while ( my $gff_line = $brs_iterator->() ) {
        next COORD unless ref $gff_line eq 'HASH';
        ++$counter;
        $lc->($gff_line) if $opt{locus_tag};

        my ($gstart, $gend, $gscore) = @{$gff_line}{qw/start end score/};
        my $glen = $gend-$gstart+1;

        for my $target (@$targets){
            my ($tstart, $tend) = @{$target}[0,1];
            my $tlen = $tend-$tstart+1;
            my $overlap = overlap($gstart, $gend, $tstart, $tend);
            $stat->add_data(
                ($opt{source} ? ($overlap / $glen) : 1) * 
                ($opt{dest}   ? ($overlap / $tlen) : 1) * 
                $gscore
            );
        }
    }

    if ($counter) {
        return {
            score => $stat->mean(),
            std   => $stat->standard_deviation(),
            var   => $stat->variance(),
            n     => $counter,
            ($opt{locus_tag} ? $lc->() : ())
        };
    }
}



##########################################################
# average_scores
# for ranges returned by brs_iterator, calculate average/stddev

sub average_scores {
    my ($brs_iterator) = @_;

    my ( $score_avg, $score_std, $score_var, $score_count ) = ( 0, 0, 0, 0 );

COORD:
    while ( my $gff_line = $brs_iterator->() ) {

        next COORD unless ref $gff_line eq 'HASH';

        my $previous_score_avg = $score_avg;

        $score_avg = $score_avg
            + ( $gff_line->{score} - $score_avg ) / ++$score_count;

        $score_std
            = $score_std
            + ( $gff_line->{score} - $previous_score_avg )
            * ( $gff_line->{score} - $score_avg );

        $score_var = $score_std / ( $score_count - 1 )
            if $score_count > 1;

    }

    if ($score_count) {
        return {
            score => $score_avg,
            std   => sqrt($score_var),
            var   => $score_var,
            n     => $score_count,
        };
    }
}

#######################################################################
# microarray

sub microarray_average{
    # print STDERR "+\n";
    my ($brs_iterator, %opt) = @_;
    my $targets = $opt{targets};

    my ( $score_avg, $score_std, $score_var, $score_count ) = ( 0, 0, 0, 0 );


    if (@$targets != 1){
        die "should only be one target";
    }

    my ($tstart, $tend) = @{$targets->[0]}[0,1];

    my @lines;
    COORD:
    while ( my $gff_line = $brs_iterator->() ) {
        next COORD unless ref $gff_line eq 'HASH';
        push @lines, $gff_line;
    }

    if (@lines == 1) {
        return {
            score => $lines[0]->{score},
            n => 1,
        }
    }
    elsif (@lines > 1){
        my @window = (0) x ($tend - $tstart + 1);
        for my $line (@lines) {
            # min and max coord within the window
            my $min = max($line->{start}, $tstart);
            my $max = min($line->{end}, $tend);

            my $relmin = $min - $tstart;
            my $relmax = $max - $tstart;

            @window[$relmin .. $relmax] = (1) x ($relmax - $relmin + 1);
            #print STDERR "$min ($relmin) $max ($relmax)\n";
        }

        my $total_coverage = sum @window;
        die Dumper \@lines, [$tstart, $tend] if $total_coverage == 0;

        my $score = (sum map {
            overlap($_->{start}, $_->{end}, $tstart, $tend) * $_->{score} 
        } @lines) / ($total_coverage);

        return {
            score => $score,
            n => scalar(@lines),
            total_coverage => $total_coverage, 
        }
    }
}

##########################################################
# binary_range_search, make_gff_iterator, read_gff 
# same as overlap_gff.pl

sub binary_range_search {
    my %options = @_;

    my $targets = $options{targets} || croak 'Need a targets parameter';
    my $queries = $options{queries} || croak 'Need a queries parameter';

    my ( $low, $high ) = ( 0, $#{$queries} );     #what does this do?
    my @iterators = ();

TARGET:
    for my $range (@$targets) {

    RANGE_CHECK:
        while ( $low <= $high ) {

            my $try = int( ( $low + $high ) / 2 );

            $low = $try + 1, next RANGE_CHECK
                if $queries->[$try]{end} < $range->[0];
            $high = $try - 1, next RANGE_CHECK
                if $queries->[$try]{start} > $range->[1];

            my ( $down, $up ) = ($try) x 2;
            my %seen = ();

            my $brs_iterator = sub {

                if (    $queries->[ $up + 1 ]{end} >= $range->[0]
                    and $queries->[ $up + 1 ]{start} <= $range->[1]
                    and !exists $seen{ $up + 1 } )
                {
                    $seen{ $up + 1 } = undef;
                    return $queries->[ ++$up ];
                }
                elsif ( $queries->[ $down - 1 ]{end} >= $range->[0]
                    and $queries->[ $down - 1 ]{start} <= $range->[1]
                    and !exists $seen{ $down - 1 }
                    and $down > 0 )
                {
                    $seen{ $down - 1 } = undef;
                    return $queries->[ --$down ];
                }
                elsif ( !exists $seen{$try} ) {
                    $seen{$try} = undef;
                    return $queries->[$try];
                }
                else {
                    return;
                }
            };
            push @iterators, $brs_iterator;
            next TARGET;
        }
    }

# In scalar context return master iterator that iterates over the list of range iterators.
# In list context returns a list of range iterators.
    return wantarray
    ? @iterators
    : sub {
        while (@iterators) {
            if ( my $gff = $iterators[0]->() ) {
                return $gff;
            }
            shift @iterators;
        }
        return;
    };
}

sub gff_read {
    return [] if $_[0] =~ m/^
                            \s*
                            \#+
                           /mx;

    my ($seqname, $source, $feature, $start, $end,
        $score,   $strand, $frame,   $attribute
    ) = split m/\t/xm, shift || return;

    $attribute =~ s/[\r\n]//mxg;

    return {
        'seqname'   => lc $seqname,
        'source'    => $source,
        'feature'   => $feature,
        'start'     => $start,
        'end'       => $end,
        'score'     => $score,
        'strand'    => $strand,
        'frame'     => $frame,
        'attribute' => $attribute
    };
}

sub make_gff_iterator {
    my %options = @_;

    my $parser     = $options{parser};
    my $file       = $options{file};
    my $GFF_HANDLE = $options{handle};

    croak
        "Need parser function reference and file name or handle to build iterator"
        unless $parser
            and ref $parser eq 'CODE'
            and (
                ( defined $file and -e $file )
                xor(defined $GFF_HANDLE and ref $GFF_HANDLE eq 'GLOB'
                        or $GFF_HANDLE eq 'ARGV'
                )
            );

    if ($file) {
        open $GFF_HANDLE, '<', $file
            or croak "Can't read $file: $!";
    }

    return sub {
        $parser->( scalar <$GFF_HANDLE> );
    };
}

=head1 NAME

 window_gff.pl - Average GFFv3 data over a sliding window or against a GFFv3 annotation file

=head1 SYNOPSIS

Run a sliding window of 50bp every 25bp interval on methylation data (attribute field assumed to be 'c=n; t=m')

 window_gff.pl --width 50 --step 25 --scoring meth --output out.gff in.gff

Average data per each locus (on score field)

 window_gff.pl --gff genes.gff --scoring average --output out.gff in.gff

Merge 'exon' features in GFF annotation file per parent ID (eg. per gene)

 window_gff.pl --gff exons.gff --scoring average --output out.gff --tag Parent --merge exon in.gff

=head1 DESCRIPTION

 This program works in two modes. It will either run a sliding window through the input GFF data,
 or load a GFF annotation file and use that to average scores in input GFF data.

 Scores calculation by default is done via fractional methylation, in which case the program
 requires that the GFF attribute field contain both 'c' and 't' tags (eg. 'c=n; t=m').
 Alternatively, the -a, --average switch will make the program average the scores in the score field.

 An additional switch, -m, --merge, will take the name of a feature in the annotation GFF file.
 It will then try to average out all the features that have a common parent locus.
 This assumes that, for example if averaging exon features: the user must use '--merge exon',
 '--tag Parent' (in the case of Arabidopsis' annotations), which will fetch a gene ID like ATG01010.1.
 It is assumed that the Parent feature of this exon is ATG010101.

=head1 OPTIONS

 window_gff.pl [OPTION]... [FILE]...

=over
 
=item -w, --width       

sliding window width (default: 50, integer)

=item -s, --step        

sliding window interval (default: 50, integer)

=item -c, --scoring     

score computation scheme.  can be: meth (default), average, sum
weighed_average, half_weighed_average, or microarray. See SCORING section for
enlightenment.

=item -m, --merge       

merge this feature as belonging to same locus. (default: no, string [eg: exon])

=item -n, --no-sort     

GFFv3 data assumed sorted by start coordinate. (default: no)

=item -k, --no-skip     

print windows or loci for which there is no coverage  (deftaul: no)

=item -g, --gff         

GFFv3 annotation file. Not compatible with --width or --step.

=item -t, --tag         

attribute field tag in _annotation_ file from which to extract locus ID (default: ID, string)

=item -u, --collect-tag 

If this attribute is given, the attribute with the field tag in _input_ file will be collected and reported in the
output. (default: [null]. Example: 'ID').

=item -f, --feature     

overwrite GFF feature field with this label           (default: no, string)

=item -b, --absolute    

organism name to fetch chromosome lengths, implies -k (default: no, string [available: arabidopsis, rice, puffer])

=item -r, --reverse     

reverse score and counts.  (not sure what this option is... ask yvonne?)

=item -o, --output      

filename to write results to (defaults to STDOUT)

=item -v, --verbose     output perl's diagnostic and warning messages

=item -q, --quiet       supress perl's diagnostic and warning messages

=item -h, --help        print this information

=item -m, --manual      print the plain old documentation page

=back

=head1 SCORING

possible values for --scoring are:

=over

=item sum:

For each window, sum the score from each overlap and report in the score column (column 6). 

=item meth:

For each window, sum the 'n', 'c', and 't' from each overlap and report the fractional methylation c/(c+t) as the score. 

=item average:

For each window, calculate the mean, standard deviation, variance for the scores of overlapping queries.

=item weighed_average:

For each window, calculate the mean, standard deviation, variance for the FULL DILUTED scores.  Diluted means that, if
the windows are arranged like so:

 Query:    |-------------------|                 Length: x, Score: n
 Window:                |--------------------|   Length: y
 Overlap:               |------|                 Length: z

Then the score contribution of the query to the window is n * (x/z) * (y/z).  This was yvonne's idea so if it doesn't
make sense, blame her.

=item half_weighed_average:

For each window, calculate the mean, standard deviation, variance for the HALF DILUTED scores.  Diluted means that, if
the windows are arranged like so:

 Query:    |-------------------|                 Length: x, Score: n
 Window:                |--------------------|   Length: y
 Overlap:               |------|                 Length: z

Then the score contribution of the query to the window is n * (x/z).

=item microarray:

For each window, if there is a single input query window mapping to it, use its score.
If there are multiple input query windows mapping, take its weighed average. For example,

  A:score 45        B:score 123
 |----------|      |-----------| input queries
          |-----------|          target window (size 50)

If A overlaps with the target by 11 bases, and B overlaps by 15, the output
score will be ((11 * 45) + (15 * 123)) / (26).  

=head1 REVISION

 Version 0.0.1

 $Rev: 410 $:
 $Author: psilva $:
 $Date: 2010-10-06 10:52:04 -0700 (Wed, 06 Oct 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/window_gff.pl $:
 $Id: window_gff.pl 410 2010-10-06 17:52:04Z psilva $:

=head1 AUTHOR

 Pedro Silva <psilva@nature.berkeley.edu/>
 Zilberman Lab <http://dzlab.pmb.berkeley.edu/>
 Plant and Microbial Biology Department
 College of Natural Resources
 University of California, Berkeley

=head1 COPYRIGHT

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

__DATA__
#arabidopsis
chr1	30432563
chr2	19705359
chr3	23470805
chr4	18585042
chr5	26992728
chrc	154478
chrm	366924
#rice
chr01	43596771
chr02	35925388
chr03	36345490
chr04	35244269
chr05	29874162
chr06	31246789
chr07	29688601
chr08	28309179
chr09	23011239
chr10	22876596
chr11	28462103
chr12	27497214
chrc	134525
chrm	490520
#puffer
1	22981688
10	13272281
11	11954808
12	12622881
13	13302670
14	10246949
15	7320470
15_random	3234215
16	9031048
17	12136232
18	11077504
19	7272499
1_random	1180980
2	21591555
20	3798727
21	5834722
21_random	3311323
2_random	2082209
3	15489435
4	9874776
5	13390619
6	7024381
7	11693588
8	10512681
9	10554956
mt	16462
