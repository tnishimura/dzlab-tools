#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Pod::Usage;
use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/lib";
use GBUtil;

my $result = GetOptions (
    "file|f=s"   => \(my $file),
    "dir|d=s"    => \(my $dir),
    "compile|c"  => \(my $compile),
    "bed|b"      => \(my $bed),
    "ctscore|ct" => \(my $ctscore),
    "parallel|p=i" => \(my $parallel = 0),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || !$file);  

gff_to_wig(
    file    => $file,
    dir     => $dir,
    compile => $compile,
    bed     => $bed,
    ctscore => $ctscore,
    parallel => $parallel,
);

=head1 gff2wig.pl 

    file|f      (required)
    dir|d       (optional, default '.')
    compile|c   (optional, default no)
    bed|b       (optional, default no)
    ctscore|ct  (optional, default no)

=cut
