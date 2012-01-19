#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
END {close STDOUT}
$| = 1;

use FindBin;
use lib "$FindBin::Bin/lib";
use AGP;

use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) unless $opt_file;

my $agp = AGP->new(file => $opt_file);

for my $chr ($agp->groups) {
    ENTRY:
    while (defined(my $entry = $agp->next($chr))){
        next ENTRY unless $entry->{type} eq 'W';
        say join "\t", 
        $chr, 
        '.', 
        '.',
        $entry->{object_start},
        $entry->{object_end},
        '.',
        $entry->{orientation},
        '.',
        "ID=" . $entry->{component_id};
    }
}

=head1 OPTIONS

=over

=item  <file> 

=for Euclid
    file.type:        readable

=back

=cut
