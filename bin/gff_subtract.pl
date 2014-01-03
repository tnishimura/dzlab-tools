#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/../lib";
use GFF::Parser;
use GFF;

END {close STDOUT}
$| = 1;

if (2 != @ARGV){
    say <<"END";
$0 - subtract Y's scores from X's and output new GFF
usage: $0 X.gff Y.gff > result.gff

NOTE: X.gff and Y.gff are assumed to be sorted with a command like
 sort -f -k1,1 -k4,4n X.gff -o X.gff
END
    exit 1;
}

my @files = @ARGV;



my $parser_left = GFF::Parser->new(file => $files[0], normalize => 0);
my $parser_right = GFF::Parser->new(file => $files[1], normalize => 0);

LEFT:
while (defined(my $left = $parser_left->next)){
    RIGHT:
    while (defined(my $right = $parser_right->next)){
        if(GFF::start_position_equal($left,$right)){
            say join "\t", $left->sequence, $left->source // '.', $left->feature // '.', 
            $left->start, $left->end, ($left->score - $right->score), $left->strand // '.', $left->frame //'.',
            $left->attribute_string // '.';
            last RIGHT;
        }
        elsif(GFF::start_position_lessthan($left,$right)){
            $parser_right->rewind($right);
            say $left->to_string;
            last RIGHT;
        }
        elsif(GFF::start_position_lessthan($right,$left)){
            say join "\t", $right->sequence, $right->source // '.', $right->feature // '.', 
            $right->start, $right->end, (0 - $right->score), $right->strand // '.', $right->frame //'.',
            $right->attribute_string;
        }
    }
}
while (defined(my $right = $parser_right->next)){
    $right->score(-$right->score);
    say $right->to_string;
}
