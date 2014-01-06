#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use List::Util qw/first max min/;
use File::Basename qw/basename/;
use Pod::Usage;
use Getopt::Long;
use Scalar::Util qw/looks_like_number/;

my $result = GetOptions (
    "help"    => \my $help,
    "use-fullpath|f" => \my $fullpath,
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || $help || !@ARGV);  

#######################################################################

my %loaded = map { $_, load_avg($_) } @ARGV;

my @indices = check_all_same(\%loaded);

say join "\t", 'bin', map { $fullpath ? $_ : basename($_) } @ARGV;

for my $i (@indices) {
    my @line = ($i);
    for my $f (@ARGV) {
        push @line, $loaded{$f}{$i};
    }
    say join "\t", @line;
}

#######################################################################

sub check_all_same{
    my $all_scores = shift;
    my %file2prop;
    while (my ($file,$sc) = each %$all_scores) {
        $file2prop{$file} = ends_properties($sc);
    }
    if (keys(%{ {map { $_, 1 } values %file2prop} }) != 1){
        say STDERR "not all ends avg files have same min/max/step";
        for my $file (sort keys %$all_scores) {
            say STDERR "$file: $file2prop{$file}";
        }
        exit 1;
    }
    return ends_indices(first {defined} (values %$all_scores));
}

sub ends_properties{
    my $scores = shift;
    my @indices = ends_indices($scores);
    my $count = scalar @indices;
    my $min = min(@indices);
    my $max = max(@indices);
    my $step = ($max - $min) / ($count - 1);
    return "$min to $max by $step";
}

sub ends_indices{
    my $scores = shift;
    return sort {$a<=>$b} keys %$scores;
}

# avg := { index => score }
sub load_avg{
    my $file = shift;
    my %scores;
    open my $fh, '<:crlf', $file;
    <$fh>;
    while (defined(my $line = <$fh>)){
        chomp $line;
        my ($pos, $sc) = split /\t/, $line;
        # if (! looks_like_number($pos) || ! looks_like_number($sc)){
        if (! looks_like_number($pos)){
            die "line $. of $file does not look like part of and ends average file: \n$line\n";
        }
        $scores{$pos} = $sc;
    }

    close $fh;
    return \%scores;
}

=head1 ends_avg_mashup.pl 

You can haz combined ends average files:

 ends_avg_mashup.pl bs-sequel-output-dir/ends/*CG*.avg > CG-ends.avg

=cut

