package Tree::Range::TestUtils;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use Carp;
use autodie;
use Tree::Range;
use List::Util qw/min max/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(sort_results gen_rand_tree search_linear randrange);

my $max = 100000;
my $N   =   1000;

sub randcoord{ return int(rand($max))-int($max/2); }
sub randlen{ return int(rand($max/100)) }

sub randrange{
    my $start = randcoord();
    return ($start, $start + randlen());
}

sub gen_rand_tree{
    my $num = shift // $N;
    my $sr = Tree::Range->new();
    for (1..$num){
        my ($start, $end) = randrange();
        my $label = "$start -> $end";
        $sr->add( $start,$end, $label);
    }
    $sr->finalize();
    return $sr;
}


sub sort_results{
    return [sort {$a->{item} cmp $b->{item}} @_];
}

1;
