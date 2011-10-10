package GFF::Statistics;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use GFF::Parser;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(getstats);

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
# 

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
#######################################################################
# statistics

sub getstats{
    my $singlec = shift;
    my %nuclear_ct;
    my %chr_ct;
    my %mit_ct;
    my $nuclear_methyl = make_averager();
    my $chr_methyl = make_averager();
    my $mit_methyl = make_averager();

    my $parser = GFF::Parser->new(file => $singlec);
    PARSE:
    while (defined(my $gff = $parser->next())){
        my $c = $gff->get_column('c');
        my $t = $gff->get_column('t');
        next PARSE if ! (defined $c && defined $t);
        my $ct = $c+$t;

        my $methyl = ($ct == 0) ? 0 : $c / ($ct);

        given ($gff->sequence){
            when (/chrc/i){
                $chr_methyl->($methyl);
                ++$chr_ct{$ct};
            }
            when (/chrm/i){
                $mit_methyl->($methyl);
                ++$mit_ct{$ct};
            }
            when (/chr\d+/i){
                $nuclear_methyl->($methyl);
                ++$nuclear_ct{$ct};
            }
            default{
                next PARSE;
            }
        }
    }

    return {
        chr_ct_median   => histmedian(\%chr_ct),
        mit_ct_median   => histmedian(\%mit_ct),
        nuc_ct_median   => histmedian(\%nuclear_ct),
        chr_ct_mean     => histmean(\%chr_ct),
        mit_ct_mean     => histmean(\%mit_ct),
        nuc_ct_mean     => histmean(\%nuclear_ct),
        chr_methyl_mean => $chr_methyl->(),
        mit_methyl_mean => $mit_methyl->(),
        nuc_methyl_mean => $nuclear_methyl->(),
        coverage        => sumhists(\%chr_ct, \%mit_ct, \%nuclear_ct),
    };
}


1;

