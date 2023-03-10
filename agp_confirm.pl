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
use FastaReader;
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

$logger->info("Slurping Assembled");
my $f = FastaReader->new(file => $opt_scaffold, slurp => 1);

#######################################################################
# Meat

$logger->info("Slurping Scaffold (the pieces)");
my $contigs = FastaReader->new(file => $opt_fasta, slurp => 1);

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

            my $contig = $contigs->get($component->{component_id});
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

