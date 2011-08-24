#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Util;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if ! $opt_list || ! $opt_input;

open my $fh, '<', $opt_list;
my %id;
while (defined(my $line = <$fh>)){
    chomp $line;
    if($line =~ /(\w+)/){
        $id{$1} = 1;
    }
}
close $fh;


open my $gff_fh, '<', $opt_input;
while (defined(my $line = <$gff_fh>)){
    my $gff = parse_gff($line);
    my $tag = $gff->get_column($opt_attribute_id);
    if (defined $tag){
        if ($opt_keep && exists $id{$tag}){
            print $line;
        }
        elsif (! $opt_keep && ! exists $id{$tag}){
            print $line;
        }
        # elsif (! $opt_keep && exists $id{$tag}){
        #     # delete it by doing nothing
        # }
    }
    elsif ($opt_keep){ # keep if no tag and default is to keep
        print $line;
    }
    # elsif (! $opt_keep){
    #     # delete it by doing nothing
    # } 
}




=head1 NAME

filter_gff_list.pl - Delete gff lines with certain attribute/locus id.

=head1 SYNOPSIS

Usage examples:

 filter_gff_list.pl -l list input.gff 

=head1 OPTIONS

=over

=item  <input> 

=for Euclid
    input.type:        readable

=item  -a <attr> | --attribute-id <attr>

=for Euclid
    attr.default:     'ID'

=item  -l <file> | --list <file>

File with list of attributes.

=for Euclid
    file.type:        readable

=item -k | --keep

Invert the operation by KEEPing the gff lines with ID's in list. (Default is to delete, not keep)

=back

=cut
