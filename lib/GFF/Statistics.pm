package GFF::Statistics;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use GFF::Parser;
use GFF::Parser::Splicer;
use Scalar::Util qw/looks_like_number/;
use List::MoreUtils qw/all/;
use Statistics::Descriptive;
use YAML qw/Bless Dump/;
use Tie::IxHash;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(gff_info methyl_stats gff_detect_width);
our @EXPORT = qw();

our @wanted_percentiles = qw/.05 .25 .50 .75 .95/;
our $DEBUG = $ENV{DZDEBUG};

#######################################################################
#                       gff_info
#######################################################################

=head2 gff_info("file1.gff", "file2.gff")

Compute basic stats about gff file, including list of features/sequences,
window length statistics.  Returns a hashref of sequences to a hashref with the
following keys:

   count mean standard_deviation min median max number_of_sequences sequences number_of_features features

=cut
sub gff_info{
    my $file = shift || croak "gff_info - arg error";

    my $count = 0;
    my $lengths = Statistics::Descriptive::Full->new();
    my %file_yaml;

    my $parser = GFF::Parser->new(file => $file, normalize => -1);
    while (defined(my $gff = $parser->next())){
        my ($seq, $feature) = ($gff->sequence(), $gff->feature());

        $file_yaml{sequences}{$seq}++;
        $file_yaml{features}{$feature}++;
        $file_yaml{count}++;

        $lengths->add_data($gff->end - $gff->start + 1);

        if ($count++ % 25000 == 0){
            say STDERR $count;
        }
    }

    $file_yaml{mean}                = sprintf "%.1f", $lengths->mean();
    $file_yaml{standard_deviation}  = sprintf "%.1f", $lengths->standard_deviation();
    $file_yaml{min}                 = $lengths->min();
    $file_yaml{median}              = $lengths->median();
    $file_yaml{max}                 = $lengths->max();
    $file_yaml{number_of_sequences} = scalar(keys %{$file_yaml{sequences}});
    $file_yaml{number_of_features}  = scalar(keys %{$file_yaml{features}});

    Bless(\%file_yaml)->keys([qw/
        count mean standard_deviation min median max number_of_sequences sequences number_of_features features
        /]);

    return \%file_yaml;
}


#######################################################################
#                       gff_detect_width
#######################################################################

=head2 gff_detect_width "file", N

Try to guess the window width of a file from its first N lines, default 100.
If uneven or not determinable, returns undef.

=cut

sub gff_detect_width{
    my $file = shift or croak "need file";
    my $numlines = shift // 100; # number of lines to read
    my $p = GFF::Parser->new(file => $file);
    my @widths;
    my $counter = 0;
    while ($counter < $numlines && defined(my $gff = $p->next())){
        $counter++;
        my ($start, $end) = ($gff->start(), $gff->end());

        croak "bad start/end column in gff_detect_width $file" 
        if (! looks_like_number($start) || !looks_like_number($end));

        my $width = $end - $start + 1;
        push @widths, $width;
    }

    if ($counter == 0){
        return;
    }
    elsif ($counter == 1){
        return $widths[0];
    }
    else{
        pop @widths; # b/c the last one can be uneven
        my $first = shift @widths;
        if (all { $first == $_ } @widths){
            return $first;
        }
        else{
            return;
        }
    }
}

#######################################################################
#                       methylation_stats
#######################################################################


#our @rownames = qw/
#nuc_ct_mean nuc_ct_median 
#chr_ct_mean chr_ct_median 
#mit_ct_mean mit_ct_median 
#nuc_methyl_mean chr_methyl_mean mit_methyl_mean 
#chr_methyl_total mit_methyl_total nuc_methyl_total
#coverage
#/;

#######################################################################
# utility

# my $c_avger = make_averager();
# $avg = $c_avger->(3.14);
sub make_averager{
    my ($current_average, $count) = @_ == 2 ? @_ : (undef, 0);
    return sub{
        my $newval = shift;
        if (defined $newval){
            if (!defined $current_average){
                ($current_average, $count) = ($newval, 1);
            }
            else{
                $current_average = ($current_average * $count + $newval) / ($count + 1);
                $count += 1;
            }
        }
        if ($current_average){
            return $current_average;
        }
        return;
    }
}

#######################################################################
# * histogram is { value => count } hashref.

# simple sum
sub sum{
    my $total = 0;
    for my $v (@_) {
        $total += $v;
    }
    return $total;
}

# convert list of values to histogram
sub tohist{
    my %accum;
    for my $val (@_) {
        $accum{$val}++;
    }
    return %accum;
}


sub histmean{
    my $hist = shift;
    if (keys %$hist == 0){ return 'na'; }

    my $total_bin_count = sum values %$hist;
    return sumhists($hist)/$total_bin_count;
}

sub histstd{
    my $hist = shift;
    if (keys %$hist == 0){ return 'na'; }

    my $total_bin_count = sum values %$hist;
    my $mean = histmean($hist);

    my $accum=0;
    while (my ($bin,$bincount) = each %$hist) {
        $accum += $bincount * ($bin - $mean) * ($bin - $mean);
    }
    return sqrt($accum/$total_bin_count);
}

sub histpercentiles{
    my $hist = shift;
    #my @wanted_percentiles = @_;
    my $num_percentiles = @wanted_percentiles;

    if (keys(%$hist) == 0){ 
        return [('na') x $num_percentiles]; 
    }

    #tie my %percentiles, "Tie::IxHash", @wanted_percentiles;
    my %percentiles;

    my $total_bin_count = sum values %$hist;
    my $bins_counted = 0;

    my @wanted_percentiles_copy = @wanted_percentiles;

    BIN:
    for my $bin (sort {$a<=>$b} keys %$hist){
        last BIN if (!@wanted_percentiles_copy);
        $bins_counted += $hist->{$bin};
        if ($bins_counted > $total_bin_count * $wanted_percentiles_copy[0]){
            $percentiles{shift @wanted_percentiles_copy} = $bin;
        }
    }
    for my $p (@wanted_percentiles) {
        if (! exists $percentiles{$p}){
            $percentiles{$p} = 'na';
        }
    }
    if (keys(%percentiles) == $num_percentiles){
        return [@percentiles{sort keys %percentiles}];
    }
    else{
        croak "BUG with histpercentiles, please report with this report: " . Dumper \%percentiles;
    }
}

sub histmedian{
    my $hist = shift;
    if (keys %$hist == 0){ return 'na'; }

    my $total_bin_count = sum values %$hist;
    my $bins_counted = 0;
    my @bins = sort {$a<=>$b} keys %$hist;

    for my $bin (@bins) {
        $bins_counted += $hist->{$bin};
        if ($bins_counted > $total_bin_count / 2){
            return $bin;
        }
    }
    return $bins[-1];
}

sub sumhists{
    my $total = 0;
    if (@_ == 0){ return 0 }
    for my $hist (@_) {
        while (my ($bin,$bin_count) = each %$hist) {
            $total += $bin * $bin_count;
        }
    }
    return $total;
}

sub combine_hists{
    my %hist;
    for my $h (@_) {
        while (my ($bin,$count) = each %$h) {
            $hist{$bin} += $count;
        }
    }
    # croak Dumper \@_, \%hist;
    return \%hist;
}

#######################################################################
# statistics

# returns multilevel hash as first value:
# { combined | cg | chg | chh }
#   { chr | nuclear | mit | total }
#     { ct_percentiles | coverage | methyl_avg | c_count | overall_methylation | t_count | line_count | ct_mean}
# and a formatted table as the second
#
# my ($stats, $table) = methyl_stats(cg => 'cg.gff', chg => 'chg.gff', chh => 'chh.gff');
sub methyl_stats{
    my %files = @_;
    # say Dumper \%files;
    my %stats;
    for my $context (sort keys %files) {
        my $file = $files{$context};
        $stats{$context} = collect_stats($file);
    }
    $stats{combined} = combine_stats(values %stats);

    my @columns = qw/ line_count methyl_avg c_count t_count coverage overall_methylation ct_mean ct_percentiles /;

    my @output;

    push @output, join "\t", 'context', 'type', @columns[0..$#columns-1], 
        map { "ct_" . $_ * 100 . "%" } @wanted_percentiles;

    for my $context (sort keys %stats) {
        for my $type (qw/nuclear chr mit total/){
            my @row = ($context, $type);

            delete $stats{$context}{$type}{ct_hist};

            for my $column (@columns) {
                my $val = $stats{$context}{$type}{$column};
                if (ref $val eq 'ARRAY'){
                    push @row, @$val;
                }
                else {
                    push @row, round_if_decimal($val);
                }
            }
            push @output, join "\t", @row;
        }
    }

    return \%stats, join "\n", @output;
}

sub round_if_decimal{
    my $num = shift;
    $num =~ /^\d+\.\d+$/ ?  sprintf("%0.4f", $num) : $num;
}

sub combine_stats{
    my @stats = @_;
    #say Dumper \@stats;
    #my @wanted_percentiles = qw/.05 .25 .50 .75 .95/;
    my @types = qw/nuclear chr mit/;
    my %combined;
    for my $type (qw/nuclear chr mit total/) {
        $combined{$type}{ct_hist} = combine_hists(map { $_->{$type}{ct_hist} }  @stats);
        $combined{$type}{c_count}    = sum(map { $_->{$type}{c_count} }  @stats);
        $combined{$type}{t_count}    = sum(map { $_->{$type}{t_count} }  @stats);
        $combined{$type}{coverage}   = sum(map { $_->{$type}{coverage} }  @stats);
        $combined{$type}{line_count} = sum(map { $_->{$type}{line_count} }  @stats);
        $combined{$type}{overall_methylation} = $combined{$type}{coverage} > 0 ?  $combined{$type}{c_count} / $combined{$type}{coverage} : 'na';

        $combined{$type}{ct_percentiles} = histpercentiles($combined{$type}{ct_hist}, @wanted_percentiles);
        $combined{$type}{ct_mean}        = histmean($combined{$type}{ct_hist});

        $combined{$type}{methyl_avg} = $combined{$type}{line_count} > 0 ? (
            sum map 
            { 
                $_->{$type}{methyl_avg} * $_->{$type}{line_count} 
            } @stats
        ) / $combined{$type}{line_count} 
        : 0;
    }
    return \%combined;
}

# return {
#    nuclear 
#    chr 
#    mit 
#    total
# }{
#    methyl_avg          # avg of scores column 
#    ct_hist 
#    c_count 
#    t_count 
#    line_count          # num lines
#    ct_percentiles      # percentile of c+t for each line 
#    ct_mean             # mean of c+t for each line 
#    coverage            # c_count + t_count (or sumhist(ct_hist))
#    overall_methylation # c_count / (c_count + t_count)
# }

sub collect_stats{
    my $singlec = shift;
    say STDERR "collect_stats($singlec)" if $DEBUG;
    my %stats;

    # my %methyl_averager = (
    #     nuclear => make_averager(),
    #     mit => make_averager(),
    #     chr => make_averager(),
    # );
    # $stats{nuclear}{methyl_avg} = undef;
    # $stats{chr}{methyl_avg}     = undef;
    # $stats{mit}{methyl_avg}     = undef;
    # $stats{nuclear}{ct_hist}    = {};
    # $stats{chr}{ct_hist}        = {};
    # $stats{mit}{ct_hist}        = {};
    # $stats{nuclear}{c_count}    = 0;
    # $stats{chr}{c_count}        = 0;
    # $stats{mit}{c_count}        = 0;
    # $stats{nuclear}{t_count}    = 0;
    # $stats{chr}{t_count}        = 0;
    # $stats{mit}{t_count}        = 0;
    # $stats{nuclear}{line_count} = 0; # redundant but more clear.
    # $stats{chr}{line_count}     = 0;
    # $stats{mit}{line_count}     = 0;

    my %methyl_averager;
    my @all_types = qw/nuclear chr mit unknown/; 
    for my $type (@all_types) {
        $methyl_averager{$type} = make_averager();
        $stats{$type}{methyl_avg} = undef;
        $stats{$type}{ct_hist}    = {};
        $stats{$type}{c_count}    = 0;
        $stats{$type}{t_count}    = 0;
        $stats{$type}{line_count} = 0; # redundant but more clear.
    }

    my $parser = GFF::Parser::Splicer->new(file => $singlec, columns => [qw/seq c t/]);

    my $counter = 0;
    PARSE:
    while (defined(my $gff = $parser->next())){
        say STDERR $counter if ++$counter % 50000 == 0 and $DEBUG;
        my ($seq, $c, $t) = @$gff;
        next PARSE if ! (defined $c && defined $t);

        my $ct = $c+$t;
        my $methyl = ($ct == 0) ? 0 : $c / ($ct);

        my $type = $seq =~ /chrc|chrpt/i ? 'chr' :
                   $seq =~ /chrm/i       ? 'mit' : 
                   $seq =~ /chr\d+/i     ? 'nuclear' : 
                   $seq =~ /unknown/i    ? 'unknown' : 
                   undef;

                   # die if ! defined $type;

        next PARSE if ! $type;

        if (! exists $methyl_averager{$type}){
            say $type;
        }
        $methyl_averager{$type}->($methyl);
        $stats{$type}{ct_hist}{$ct}++;
        $stats{$type}{c_count} += $c;
        $stats{$type}{t_count} += $t;
        $stats{$type}{line_count}++;
    }

    for my $type (@all_types) {
        $stats{$type}{methyl_avg} = $methyl_averager{$type}->() // 0;
    }
    # $stats{nuclear}{methyl_avg} = $methyl_averager{nuclear}->() // 0;
    # $stats{chr}{methyl_avg}     = $methyl_averager{chr}->() // 0;
    # $stats{mit}{methyl_avg}     = $methyl_averager{mit}->() // 0;

    #######################################################################
    # total - never includes unknown
    $stats{total}{line_count} = $stats{nuclear}{line_count} + $stats{mit}{line_count} + $stats{chr}{line_count};

    $stats{total}{methyl_avg} = $stats{total}{line_count} == 0 ? 0 : (
        $stats{nuclear}{methyl_avg} * $stats{nuclear}{line_count} + 
        $stats{mit}{methyl_avg} * $stats{mit}{line_count} + 
        $stats{chr}{methyl_avg} * $stats{chr}{line_count}
    ) / $stats{total}{line_count};

    $stats{total}{ct_hist} = combine_hists(
        $stats{nuclear}{ct_hist},
        $stats{mit}{ct_hist},
        $stats{chr}{ct_hist},
    );
    $stats{total}{c_count} = $stats{nuclear}{c_count} + $stats{mit}{c_count} + $stats{chr}{c_count};
    $stats{total}{t_count} = $stats{nuclear}{t_count} + $stats{mit}{t_count} + $stats{chr}{t_count};

    # for my $type (exists $stats{unknown} ? qw/chr mit nuclear unknown total/ : qw/chr mit nuclear total/){
    for my $type (qw/chr mit nuclear total/){
        $stats{$type}{ct_percentiles} = histpercentiles($stats{$type}{ct_hist}, @wanted_percentiles);
        $stats{$type}{ct_mean}        = histmean($stats{$type}{ct_hist});
        $stats{$type}{coverage}       = $stats{$type}{c_count} + $stats{$type}{t_count};
        $stats{$type}{overall_methylation} = 
          ($stats{$type}{coverage} > 0) ?  ($stats{$type}{c_count} / $stats{$type}{coverage}) : 'na';
    }

    return \%stats;
}

1;

