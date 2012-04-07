package GFF::Statistics;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use GFF::Parser;
use Scalar::Util qw/looks_like_number/;
use List::MoreUtils qw/all/;
use Statistics::Descriptive;
use YAML qw/Bless Dump/;
use Tie::IxHash;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(gff_info methylation_stats gff_detect_width);
our @EXPORT = qw();

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
        return 0;
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
        return 0;
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

sub make_averager{
    my ($current_average, $count) = (undef, 0);
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
        return $current_average // 'na';
    }
}

#######################################################################

sub sum{
    my $total = 0;
    for my $v (@_) {
        $total += $v;
    }
    return $total;
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
    if (keys %$hist == 0){ return 'na'; }

    my @wanted_percentiles = @_;

    tie my %percentiles, "Tie::IxHash", @wanted_percentiles;

    my $total_bin_count = sum values %$hist;
    my $bins_counted = 0;

    BIN:
    for my $bin (sort {$a<=>$b} keys %$hist){
        last BIN if (!@wanted_percentiles);
        $bins_counted += $hist->{$bin};
        if ($bins_counted > $total_bin_count * $wanted_percentiles[0]){
            $percentiles{shift @wanted_percentiles} = $bin;
        }
    }
    return \%percentiles;
    if (keys %percentiles == @wanted_percentiles){
        return \%percentiles;
    }
    else{
        return 'na';
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

sub tohist{
    my %accum;
    for my $val (@_) {
        $accum{$val}++;
    }
    return %accum;
}

# histogram inter-quartile range (from 25% to 75%)
sub hist_iqr{
    my ($hist, $percentile_aref) = @_;
    return {} if keys %$hist == 0;

    my ($percentile25, $percentile75) = @{$percentile_aref}[1,3];
    return {
        map {
            $_ => $hist->{$_}
        }
        grep {
            $percentile25 <= $_ && $_ <= $percentile75
        }
        keys %$hist
    };
}
#######################################################################
# statistics

sub methylation_stats{
    my $singlec = shift;
    my @wanted_percentiles = @_;
    if (@wanted_percentiles == 0){
        @wanted_percentiles = qw/.05 .25 .50 .75 .95/;
    }

    my %nuclear_ct;
    my %chr_ct;
    my %mit_ct;
    my $nuclear_methyl = make_averager();
    my $chr_methyl = make_averager();
    my $mit_methyl = make_averager();
    my ($mit_c, $mit_t, $chr_c, $chr_t, $nuc_c, $nuc_t) = (0) x 6;

    my $parser = GFF::Parser->new(file => $singlec);

    my $counter = 0;
    PARSE:
    while (defined(my $gff = $parser->next())){
        say STDERR $counter if $counter++ % 50000 == 0;
        my $c = $gff->get_column('c');
        my $t = $gff->get_column('t');
        next PARSE if ! (defined $c && defined $t);
        my $ct = $c+$t;

        my $methyl = ($ct == 0) ? 0 : $c / ($ct);

        given ($gff->sequence){
            when (/chrc/i){
                $chr_methyl->($methyl);
                ++$chr_ct{$ct};
                $chr_c += $c;
                $chr_t += $t;
            }
            when (/chrm/i){
                $mit_methyl->($methyl);
                ++$mit_ct{$ct};
                $mit_c += $c;
                $mit_t += $t;
            }
            when (/chr\d+/i){
                $nuclear_methyl->($methyl);
                ++$nuclear_ct{$ct};
                $nuc_c += $c;
                $nuc_t += $t;
            }
            default{
                next PARSE;
            }
        }
    }

    # for my $bin (sort { $a <=> $b } keys %nuclear_ct){
    #     say "$bin => $nuclear_ct{$bin}";
    # }

    my $chr_ct_percentiles = histpercentiles(\%chr_ct, @wanted_percentiles);
    my $mit_ct_percentiles = histpercentiles(\%mit_ct, @wanted_percentiles);
    my $nuc_ct_percentiles = histpercentiles(\%nuclear_ct, @wanted_percentiles);
    #my %chr_ct_iqr = hist_iqr(\%chr_ct, $chr_ct_percentiles);
    #my %mit_ct_iqr = hist_iqr(\%mit_ct, $mit_ct_percentiles);
    #my %nuc_ct_iqr = hist_iqr(\%nuclear_ct, $nuc_ct_percentiles);

    return {
        file => $singlec,
        nuc_c => $nuc_c, nuc_t => $nuc_t,
        mit_c => $mit_c, mit_t => $mit_t,
        chr_c => $chr_c, chr_t => $chr_t,

        chr_ct_percentiles => $chr_ct_percentiles, 
        mit_ct_percentiles => $mit_ct_percentiles,
        nuc_ct_percentiles => $nuc_ct_percentiles,

        chr_ct_mean      => histmean(\%chr_ct),
        mit_ct_mean      => histmean(\%mit_ct),
        nuc_ct_mean      => histmean(\%nuclear_ct),

        #chr_ct_std       => histstd(\%chr_ct),
        #mit_ct_std       => histstd(\%mit_ct),
        #nuc_ct_std       => histstd(\%nuclear_ct),

        chr_methyl_mean  => $chr_methyl->(),
        mit_methyl_mean  => $mit_methyl->(),
        nuc_methyl_mean  => $nuclear_methyl->(),

        chr_methyl_total => ($chr_c + $chr_t > 0) ? ($chr_c / ($chr_c+$chr_t)) : 'na',
        mit_methyl_total => ($mit_c + $mit_t > 0) ? ($mit_c / ($mit_c+$mit_t)) : 'na',
        nuc_methyl_total => ($nuc_c + $nuc_t > 0) ? ($nuc_c / ($nuc_c+$nuc_t)) : 'na',

        methyl_total     => ($chr_c + $chr_t + $mit_c + $mit_t + $nuc_c + $nuc_t > 0) ?  ($chr_c + $mit_c + $nuc_c)/($chr_c + $chr_t + $mit_c + $mit_t + $nuc_c + $nuc_t) : 'na',
            
        coverage         => sumhists(\%chr_ct, \%mit_ct, \%nuclear_ct),
    }, 
    \%nuclear_ct,
    \%chr_ct,
    \%mit_ct,
    ;
}

1;

