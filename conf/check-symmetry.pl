#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
END {close STDOUT}
$| = 1;
use Test::More qw(no_plan);
use Test::Exception;

use YAML qw/LoadFile/;

my ($forward, $backward) = @ARGV;
if (! $forward || ! $backward){
    die "usage: $0 forward.alignment backward.alignment";
}
 
my %back = %{LoadFile($forward)};
my %forward = %{LoadFile($backward)};

my @back_chr = sort keys %{$back{alignment}};
my @forward_chr = sort keys %{$forward{alignment}};

is_deeply(\@back_chr, \@forward_chr, "same chromosomes");

for my $chr (@forward_chr) {
    my $numfrags = scalar(@{$forward{alignment}{$chr}});
    is($numfrags, scalar(@{$back{alignment}{$chr}}), "same number of $chr fragments");

    for my $i (0 .. $numfrags - 1) {
        is(
            $forward{alignment}{$chr}[$i][0], 
            $back{alignment}{$chr}[$i][2], 
            "$chr frag $i $forward{alignment}{$chr}[$i][0] vs $back{alignment}{$chr}[$i][2]",
        );
        is(
            $forward{alignment}{$chr}[$i][1], 
            $back{alignment}{$chr}[$i][3], 
            "$chr frag $i $forward{alignment}{$chr}[$i][1] vs $back{alignment}{$chr}[$i][3]",
        );
    }
}

__DATA__

==> at-tair10-to-tair8.alignment <==
---
alignment:
  Chr1:
    -
      - 1
      - 838264
      - 1
      - 838264
    -
      - 838266

==> at-tair8-to-tair10.alignment <==
---
alignment:
  Chr1:
    -
      - 1
      - 838264
      - 1
      - 838264
    -
      - 838265
