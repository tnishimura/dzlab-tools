#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/lib";
use FastaOO;
use Fasta;
use AGP;
use DZUtil qw/reverse_complement/;

use Log::Log4perl qw/:easy/;
Log::Log4perl->easy_init({ level => $DEBUG, layout => '%d{HH:mm:ss} %.1p > %m%n' });
my $logger = get_logger();

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) if $opt_help;

#######################################################################
# AGP

$logger->info("Reading in AGP file");
my $agp = AGP->new(file => $opt_agp);

#######################################################################
# Scaffold

$logger->info("FastaOO");
my $f = FastaOO->new(file => $opt_scaffold);
#$logger->debug(Dumper $f->counts);

#######################################################################
# Meat

$logger->info("Slurping Fasta");
my $contigs = slurp_fasta($opt_fasta);

for my $group ($agp->groups()) {
    $logger->info($group);

    COMP:
    while (defined(my $component = $agp->next($group))){
        $logger->trace($component->{component_id});

        my $scaffold_chunk = $f->get($group, $component->{object_start}, $component->{object_end}, base => 1);

        if ($component->{type} eq 'N'){
            if ($scaffold_chunk !~ /^N+$/){
                $logger->logdie(Dumper $component);
            }
        } else{
            #$logger->logdie("can't handle minus strand yet") if $component->{orientation} eq '-';
            my $contig = $contigs->{$component->{component_id}};
            my $contig_substring = substr $contig, $component->{component_start} - 1, $component->{component_length};
            $contig_substring = reverse_complement($contig_substring) if $component->{orientation} eq '-';

            if (length $contig_substring != $component->{component_length} ||
                length $scaffold_chunk != $component->{component_length} || 
                $contig_substring ne $scaffold_chunk){
                $logger->logdie(Dumper $component);
            }
            else {
                $logger->trace($component->{component_id} . " was ok");
            }
        }
    }
}
$logger->info("OK");


=head1 NAME
 
agp_confirm.pl
 
=head1 SYNOPSIS

 agp_confirm.pl

=head1 OPTIONS

=over

=item  -a <file> | --agp <file>

=for Euclid
    file.type:        readable

=item  -f <file> | --fasta <file>

=for Euclid
    file.type:        readable

=item  -s <file> | --scaffold <file>

=for Euclid
    file.type:        readable

=item -h | --help

=back

=cut

