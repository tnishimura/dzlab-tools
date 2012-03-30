package Tree::Range::Overlap_PP;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(_overlap);

sub _overlap{
    my ($start1, $end1, $start2, $end2) = @_;

    if ($end1 >= $start2 && $end2 >= $start1){
        return ($end1 < $end2 ? $end1 : $end2) - ($start1 > $start2 ? $start1 : $start2) + 1;
    }
    else {
        return 0;
    }
}

1;
