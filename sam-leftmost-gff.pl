#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;

use FindBin;
use lib "$FindBin::Bin/lib";
use Sam::Parser;

use Pod::Usage;
my $file = shift or pod2usage(-verbose => 2, -noperldoc => 1);

my $p = Sam::Parser->new(file => $file);

while (defined(my $sam = $p->next)){
    next if ! $sam->mapped;
    say join("\t",
        $sam->seqid,
        ".",
        ".",
        $sam->original_leftmost,
        $sam->original_leftmost,
        ".",
        ".",
        ".",
        "ReadID=" . $sam->readid
    );
}

=head1 sam-leftmost-gff.pl 

Create gff file from left-most coordinates of sam file.

 sam-leftmost-gff.pl input.sam > output.gff

=cut

