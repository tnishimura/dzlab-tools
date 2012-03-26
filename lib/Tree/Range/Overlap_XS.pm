package Tree::Range::Overlap_XS;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Inline 'C';

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(_overlap);

1;

__DATA__
__C__
    int _overlap(int s1,int e1,int s2, int e2) {
        if (e1 >= s2 && e2 >= s1){
            return (e1 < e2 ? e1 : e2) - (s1 > s2 ? s1 : s2) + 1;
        }
        else {
            return 0;
        }
    }
