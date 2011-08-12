#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use List::Util qw/sum/;

END {close STDOUT}
$| = 1;


if (!$opt_kernel){
    pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) ;
    exit (1);
}

my @kernel = split /,/, $opt_kernel;
if (@kernel % 2 != 1){
    pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) ;
    exit(1);
}

my $norm = sum @kernel;
@kernel = map { $_ / $norm } @kernel;

my $half_width = int(@kernel / 2);
open my $fh, '<', $opt_input;

{
    my @buffer;
    
    for (1 .. $half_width) {
        push @buffer, undef;
    }
    for (1 .. $half_width) {
        my $line = scalar <$fh>;
        chomp $line if defined $line;
        push @buffer, $line;
    }
    # @buffer = (undef, undef, line1, line2)

    sub next_window{
        my $line = scalar <$fh>;
        chomp $line if defined $line;
        push @buffer, $line;
        while (@buffer > @kernel){
            shift @buffer;
        }
        return \@buffer;
    }
}

sub average{
    my @s = map { defined $_ ? (split /\t/, $_)[5] : 0} @_;
    
    return sprintf("%.4f",sum map { $s[$_] * $kernel[$_] } (0 .. $#kernel));
}

while (defined(my $window = next_window())){
    my $current = $window->[$half_width];
    last if ! defined $current;

    my @split = split /\t/, $current;
    $split[5] = average(@$window);
    say join "\t", @split;
}


=head1 NAME

gff_smooth.pl - smooth gff file by a specifiable kernel

=head1 SYNOPSIS

Usage examples:

 moving_average.pl -k ".1,.2,.4,.2,.1" input.gff > output.gff

=head1 REQUIRED ARGUMENTS

=over

=back

=head1 OPTIONS

=over

=item  -k <kernel_specs> | --kernel <kernel_specs>

comma separate list of kernel. 

=item <input>

=for Euclid
    input.type:        readable

=back

=cut
