#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use autodie;

my $output      = '-';
my $scores      = 100;
my $bin_width   = 1;

#$bin_width = $config{'bin-width'} // $bin_width;
#if (exists $config{distance}){
#    $scores = $config{distance} * 2 / $bin_width;
#}

# Grabs and parses command line options
my $result = GetOptions (
    'scores|s=i'    => \$scores,
    'bin-width|w=i' => \$bin_width,
    'output|o=s'    => \$output,
    'verbose|v'     => sub { use diagnostics; },
    'quiet|q'       => sub { no warnings; },
    'help|h'        => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'      => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result;

if ($output ne '-') {
    open my $USER_OUT, '>', $output;
    select $USER_OUT;
}

my @stats = map { []  } (1 .. $scores);

LOCUS:
while (<ARGV>) {
    $_ =~ tr/\n\r//d;
    my ($locus, @scores) = split /\t/;

    for my $k (0 .. $scores - 1){
        next if $k >= @scores;
        my $s = $scores[$k];
        next if ($s eq 'na');
        push @{$stats[$k]}, $s;
        #$stats[$k]->add_data($s);
    }

}

print join ("\t", qw/bin mean std var ste numscores 25% 50% 75%/), "\n";

for my $k (0 .. $scores - 1){
    my $vals = $stats[$k];
    my $count = @$vals;
    my $index = $k * $bin_width - int ($scores/2) * $bin_width;

    if ($count){
        my $mean = smean(@$vals);
        my $var = svar(@$vals);
        my $std = $var eq 'na' ? 'na' : sqrt($var);
        my $ste = $var eq 'na' ? 'na' : $std/sqrt($count);
        my @quartiles = quartiles(@$vals);

        printf("%d\t" . ("%s\t" x 4) . "%d\t" . ("%s\t" x 3) . "\n",
            $index, $mean, $std, $var, $ste, $count, @quartiles,
        );
    }
    else{
        printf("%d\t" . ("%s\t" x 4) . "%d\t" . ("%s\t" x 3) . "\n",
            $index, 'na', 'na', 'na', 'na', $count, 'na', 'na', 'na',
        );
    }
}


#######################################################################
# math

sub restrict_range{
    my ($val, $min, $max) = @_;
    return $val < $min ? $min : 
           $val > $max ? $max :
           $val;
}

sub quartiles{
    my @data = sort { $a<=>$b } @_;
    return @data[
    restrict_range(int(@data/4),0,$#data), 
    restrict_range(int(@data/2),0,$#data), 
    restrict_range(int(@data*3/4),0,$#data),
    ];
}

sub my_sum{
    my $total = 0;
    for (@_){$total += $_;}
    return $total;
}

sub smean{
    return my_sum(@_) / scalar(@_);
}

sub svar{
    if (@_ < 2){
        return 'na';
    }
    my $n = scalar @_;
    my $mean = smean(@_);
    return (my_sum(map { $_*$_ } @_) - $n * $mean * $mean ) / ($n - 1);
}

=head1 NAME

 average_ends.pl - Average ends analysis bins wit statistics

=head1 SYNOPSIS

 # Providing --bin-width, -w lets the program output proper, absolute bin indices
 average_ends.pl --bin-width 100 --scores 100 input.ends -o output.dat

=head1 DESCRIPTION

 Expects an ends analysis file in which each record represents one locus, followed by a set number of scores per bin.
 A bin is a portion of the locus, for example from its 5' ends to its 3' end in the case of genes.

=head1 OPTIONS

 average_ends.pl [OPTION]... [FILE]...

 -s, --scores      Number of scores per locus record
 -w, --bin-width   Width of each scores bin for calculating absolute bin indices
 -o, --output      filename to write results to (defaults to STDOUT)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page

=head1 REVISION

 Version 0.0.1

 $Rev: 249 $:
 $Author: psilva $:
 $Date: 2010-01-11 21:24:34 -0800 (Mon, 11 Jan 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/average_ends.pl $:
 $Id: average_ends.pl 249 2010-01-12 05:24:34Z psilva $:

=head1 AUTHOR

 Pedro Silva <pedros@berkeley.edu/>
 Zilberman Lab <http://dzlab.pmb.berkeley.edu/>
 Plant and Microbial Biology Department
 College of Natural Resources
 University of California, Berkeley

=head1 COPYRIGHT

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
