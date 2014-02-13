#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;
use autodie;

use Pod::Usage;
use Getopt::Long;
use Test::Deep::NoTest;

use FindBin;
use lib "$FindBin::Bin/lib";
use Sam::Parser;

my $result = GetOptions (
    "common|c=s" => \(my $output_file_common),
);
pod2usage(-verbose => 2, -noperldoc => 1) if ( ! $result );  

my ($input_file_a, $input_file_b) = @ARGV;
pod2usage(-verbose => 2, -noperldoc => 1) if ( ! $input_file_a || ! $input_file_b );  

my $output_file_a = "$input_file_a.filtered";
my $output_file_b = "$input_file_b.filtered";

my $ain = Sam::Parser->new(file => $input_file_a, skip_unmapped => 0);
my $bin = Sam::Parser->new(file => $input_file_b, skip_unmapped => 0);

open my $aout, '>', $output_file_a;
open my $bout, '>', $output_file_b;

# print headers
for my $h (@{$ain->header_lines}) { $aout->print("$h\n"); }
for my $h (@{$bin->header_lines}) { $bout->print("$h\n"); }

# my $cout;
# if ($output_file_common){
#     open $cout, '>', $output_file_common;
# }

CMP:
while (
    defined (my $a_sam = $ain->next) and 
    defined (my $b_sam = $bin->next)
) {
    if ($a_sam->readid ne $b_sam->readid) {
        die "$a_sam\n$b_sam\ndo not match";
    }

    my $a_mapped = $a_sam->mapped;
    my $b_mapped = $b_sam->mapped;
    my @a_snps = $a_mapped ? @{$a_sam->snps} : ();
    my @b_snps = $b_mapped ? @{$b_sam->snps} : ();
    # my $a_num_snps = scalar @a_snps;
    # my $b_num_snps = scalar @b_snps;
    my $a_num_snps = $a_sam->mismatches_xm;
    my $b_num_snps = $b_sam->mismatches_xm;

    if ($a_mapped && $b_mapped){
        if ($a_num_snps == 0 and $b_num_snps == 0){
            # neither reported snps, nothing of interest here.
            next;
        }
        elsif ($a_num_snps < $b_num_snps){
            $aout->print("$a_sam\n");
        }
        elsif ($a_num_snps > $b_num_snps){
            $bout->print("$b_sam\n");
        }
        # same number of snps, need more info
        elsif ($a_sam->mapq > $b_sam->mapq) {
            $aout->print("$a_sam\n");
        }
        elsif ($a_sam->mapq < $b_sam->mapq) {
            $bout->print("$b_sam\n");
        }
        # diff snps... not sure about this.\    
        elsif (! eq_deeply(\@a_snps, \@b_snps)){
            $aout->print("$a_sam\n");
            $bout->print("$b_sam\n");
        }
        else{
            next;
        }
    }
    elsif ($a_mapped){
        # if ($a_num_snps > 0){
            $aout->print("$a_sam\n");
            # }
    }
    elsif ($b_mapped){
        # if ($b_num_snps > 0){
            $bout->print("$b_sam\n");
            # }
    }
    else {
        # both unmapped
        next;
    }
}

close $aout; 
close $bout;

=head2 $sam->readid 
=head2 $sam->mapped
=head2 $sam->reverse
=head2 $sam->failed_qc

__END__

=head1 NAME

 split_on_mismatches_sam.pl [--common common.sam] left.sam right.sam

=cut

