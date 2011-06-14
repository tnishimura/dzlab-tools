#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use File::Temp;
use FindBin;
use lib "$FindBin::Bin/lib";
use Launch;

use Log::Log4perl qw/:easy/;
Log::Log4perl->easy_init({ level => $DEBUG, layout => '%d{HH:mm:ss} %.1p > %m%n' });
my $logger = get_logger();

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if (!$opt_output || ! @opt_file);

my $tmpdb = mktemp($opt_output . ".tmp.XXXXX");

launch("perl -S single_c_mashup_create.pl -d $tmpdb");

while (my $file = shift @opt_file){
    if (-f $file){
        $logger->info("adding $file");
        launch("perl -S single_c_mashup_add.pl -d $tmpdb -g $file");
    }
    else {
        my $nick = $file;
        $file = shift @opt_file;
        $logger->info("adding $file");
        launch("perl -S single_c_mashup_add.pl -d $tmpdb -g $file -n $nick");
    }
}

launch("perl -S single_c_mashup_output.pl -d $tmpdb -o $opt_output");

unlink $tmpdb;

=head1 NAME
 
 single_c_mashup.pl -o output.txt [[nickname1] file1.gff] [[nickname2] file2.gff] ...
 
=head1 SYNOPSIS

 single_c_mashup.pl -o output.txt mutant1 mutant1.single-c-CG.gff.merged mutant2 mutant2.single-c-CG.gff.merged 

=head1 OPTIONS

=over

=item -o <file> | --output <file>

=item <file>...

=back

=cut

