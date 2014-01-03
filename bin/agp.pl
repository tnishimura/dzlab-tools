#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use List::Util qw/max/;
use lib "$FindBin::Bin/lib";
#use Fasta qw/slurp_fasta format_fasta/;
use FastaReader;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use AGP;
use DZUtil qw/reverse_complement/;
use Log::Log4perl qw/:easy/;
Log::Log4perl->easy_init({ level => $DEBUG, layout => '%d{HH:mm:ss} %.1p > %m%n' });
my $logger = get_logger();

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) if $opt_help;

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh; 
}

#######################################################################
# AGP

$logger->info("Reading in AGP file");
my $agp = AGP->new(file => $opt_agp);

#######################################################################
# slurp fasta

$logger->info("Slurping Fasta");
my $contigs = FastaReader->new(file => $opt_fasta, slurp => 1);
#my $contigs = slurp_fasta($opt_fasta);

#$logger->debug("contigs: \n" . join "\n", keys %$contigs);

#######################################################################
# Meat

$logger->info("groups: " . join ",", $agp->groups());

for my $group ($agp->groups) {
    my $len = $agp->length($group);
    $logger->debug("Group $group, length $len");

    my $sequence = 'N' x $len;

    COMP:
    while (defined(my $component = $agp->next($group))){
        next COMP if $component->{type} eq 'N';
        #$logger->logdie("can't handle minus strand yet") if $component->{orientation} eq '-';
        $logger->trace($component->{component_id});

        #my $contig = $contigs->{$component->{component_id}};
        my $contig = $contigs->get($component->{component_id});
        my $contig_substring = substr $contig, $component->{component_start} - 1, $component->{component_length};

        $contig_substring = reverse_complement($contig_substring) if $component->{orientation} eq '-';
        
        substr $sequence, $component->{object_start}-1, $component->{object_length}, $contig_substring;
    }

    print FastaReader::format_fasta($group, $sequence);
}

=head1 OPTIONS

=over

=item  -a <file> | --agp <file>

=for Euclid
    file.type:        readable

=item  -f <file> | --fasta <file>

=for Euclid
    file.type:        readable

=item -o <file> | --output <file>

=for Euclid
    file.default:        '-'

=item -h | --help

=back

=cut

