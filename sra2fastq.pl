#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Path::Class;

use Pod::Usage;
@ARGV > 0 or pod2usage(-verbose => 2, -noperldoc => 1);

for my $f (@ARGV) {
    my $parent = file($f)->parent();
    system('fastq-dump', '-O', $parent, $f);
}
=head1 sra2fastq.pl 

 sra2fastq.pl file1.sra file2.sra

=cut

