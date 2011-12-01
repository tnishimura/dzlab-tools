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

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(gff_detect_width);
our @EXPORT = qw(getstats);

our @rownames = qw/
nuc_ct_mean nuc_ct_median 
chr_ct_mean chr_ct_median 
mit_ct_mean mit_ct_median 
nuc_methyl_mean chr_methyl_mean mit_methyl_mean 
chr_methyl_total mit_methyl_total nuc_methyl_total
coverage
/;

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

    return {
        chr_ct_median    => histmedian(\%chr_ct),
        mit_ct_median    => histmedian(\%mit_ct),
        nuc_ct_median    => histmedian(\%nuclear_ct),
        chr_ct_mean      => histmean(\%chr_ct),
        mit_ct_mean      => histmean(\%mit_ct),
        nuc_ct_mean      => histmean(\%nuclear_ct),
        chr_methyl_mean  => $chr_methyl->(),
        mit_methyl_mean  => $mit_methyl->(),
        nuc_methyl_mean  => $nuclear_methyl->(),
        chr_methyl_total => ($chr_c + $chr_t > 0) ? ($chr_c / ($chr_c+$chr_t)) : 'na',
        mit_methyl_total => ($mit_c + $mit_t > 0) ? ($mit_c / ($mit_c+$mit_t)) : 'na',
        nuc_methyl_total => ($nuc_c + $nuc_t > 0) ? ($nuc_c / ($nuc_c+$nuc_t)) : 'na',
        coverage         => sumhists(\%chr_ct, \%mit_ct, \%nuclear_ct),
    };
}


1;

