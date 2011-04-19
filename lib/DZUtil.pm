package DZUtil;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;
use File::Spec::Functions;
use File::Basename;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(overlap chext split_names);
our @EXPORT = qw();


=head2 chext("/etc/test.txt", "newext")

return filename with new extension

=cut

sub chext{
    my ($fullpath, $ext) = @_;
    my ($filename, $path) = fileparse($fullpath);
    $filename =~ s/\.[^.]+$/.$ext/;
    return catfile($path,$filename);
}

=head2 split_names("/home/test/path.txt", qw/chr1 chr2 chr3/)

Yields (/home/test/path-chr1.txt /home/test/path-chr2.txt /home/test/path-chr2.txt)

=cut

sub split_names{
    my ($fullpath, @splits) = @_;
    my ($filename, $path) = fileparse($fullpath);
    my $ext = ($filename =~ s/\.([^.]+)$//) ? $1 : "";

    return map { catfile($path, $filename) . "-$_.$ext" } @splits;
}

=head2 overlap([start1,end1], [start2,end2])

Return overlap between two ranges, or 0 if not overlapping.

=cut

use List::Util qw/max min/;

sub overlap{
    my ($x,$y) = @_;


    if ($x->[1] >= $y->[0] && $y->[1] >= $x->[0]){
        return min($x->[1], $y->[1]) - max($x->[0], $y->[0])  + 1;
    }
    else {
        return 0;
    }
}




1;

