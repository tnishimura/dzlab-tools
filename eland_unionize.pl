#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;

use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) if $opt_help;

if ($opt_output){
    open my $fh, '>', $opt_output;
    select $fh;
}

sub slurp_eland{
    my $file = shift;
    my %accum;
    open my $fh, '<', $file;
    while (defined(my $line = <$fh>)){
        chomp $line;
        my ($readid) = (split /\t/,$line)[0];
        $accum{$readid} = $line;
    }
    close $fh;
    return \%accum;
}


my $left = slurp_eland($opt_left);
my $right = slurp_eland($opt_right);

my %common;
for (keys %$left, keys %$right) {$common{$_}=1;}
my @common_reads = keys %common;


for my $read (@common_reads){
    my $left_read = $left->{$read};
    my $right_read = $right->{$read};
    if (! defined $left_read && ! defined $right_read){
        die "$read is not in left or right?";
    } elsif (defined $left_read && defined $right_read){
        say $left_read;
    } elsif (defined $left_read){
        say $left_read;
    } elsif (defined $right_read){
        say $right_read;
    } else{
        die "why am I here?";
    }
}


=head1 OPTIONS

=over

=item  -l <eland> | --left <eland>

=for Euclid
    eland.type: readable
    
=item  -r <eland> | --right <eland>

=for Euclid
    eland.type: readable

=item  -o <eland> | --output <eland>


=item --help | -h

=back

=cut

