#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use File::Basename;
use File::Spec::Functions;
pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) if !$opt_outdir;



if (! -d $opt_outdir){
    mkdir $opt_outdir;
}

for my $script (glob("*.pl doc/*.pod")) {
    my $out = 
    $script =~ /index\.pod$/ ?  catfile($opt_outdir, "index.html") : catfile($opt_outdir, basename($script)) . ".html";
    system("podtree2html --hr 3 --notoc $script $out");
}

=head1 OPTIONS

=over

=item  -o <dir> | --outdir <dir>

=item --help | -h

=back

=cut
