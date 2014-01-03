#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;
use List::Util qw/min max/;
use Pod::Usage;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;

=head1 gff-slicer.pl 

 gff-slicer.pl -5 -f -500 -to 500 [-s] input.gff

=cut


my $result = GetOptions (
    "five-prime|5"  => \(my $five_prime),
    "three-prime|3" => \(my $three_prime),
    "from|f=i"      => \(my $from),
    "to|t=i"        => \(my $to),
    "suffix|s"      => \(my $suffix),
);
pod2usage(-verbose => 2, -noperldoc => 1) 
if (!$result || ! defined $from || ! defined $to || ! ($five_prime xor $three_prime));  

my $parser = GFF::Parser->new(file => \*ARGV, normalize => 0);

while (defined(my $gff = $parser->next)){
    my $rev = $gff->is_reverse;
    if ($suffix){
        if (my $id = $gff->get_column("ID")){
            if ($five_prime){
                $gff->attribute_string("ID=${id}-5p_${from}_${to}");
            }
            else{
                $gff->attribute_string("ID=${id}-3p_${from}_${to}");
            }
        }
    }
    my $orig_start = $gff->start;
    my $orig_end = $gff->end;

    # |------->
    # ^
    if ($five_prime && ! $rev){ 
        $gff->start($orig_start + $from);
        $gff->end  (min($orig_end, $orig_start + $to));
    }
    # <-------|
    #         ^       
    elsif ($five_prime && $rev){
        $gff->start(max($orig_start, $orig_end - $to));
        $gff->end  ($orig_end - $from);
    }
    # |------->
    #         ^
    elsif ($three_prime && ! $rev){ 
        $gff->start(max($orig_start, $orig_end + $from));
        $gff->end  ($orig_end + $to);
    }
    # <-------|
    # ^       
    elsif ($three_prime && $rev){
        $gff->start($orig_start - $to);
        $gff->end  (min($orig_end, $orig_start - $from));
    }
    else{
        die "meow";
    }
    say $gff;
}

