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
use Config::General qw(ParseConfig);
use POSIX qw/strftime/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(read_conf fastq_read_length timestamp datestamp overlap chext split_names base_match);
our @EXPORT = qw();

=head2 chext("/etc/test.txt", "newext")

return filename with new extension

=cut

sub chext{
    my ($fullpath, $ext) = @_;
    my ($filename, $path) = fileparse($fullpath);
    $filename =~ s/\.[^.]+$//;
    return catfile($path,$filename) . ".$ext";
}

=head2 split_names("/home/test/path.txt", qw/chr1 chr2 chr3/)

Yields (/home/test/path-chr1.txt /home/test/path-chr2.txt /home/test/path-chr2.txt)

=cut

sub split_names{
    my ($fullpath, @splits) = @_;
    my ($filename, $path) = fileparse($fullpath);
    my $ext = ($filename =~ s/\.([^.]+)$//) ? $1 : "";

    return map { catfile($path, $filename) . "-$_" . ($ext ? ".$ext" : "") } @splits;
}

=head2 overlap([start1,end1], [start2,end2])

Return overlap between two ranges, or 0 if not overlapping.

=cut

use List::Util qw/max min first/;

sub overlap{
    my ($x,$y) = @_;
    my $start1 = min($x->[0], $x->[1]);
    my $end1   = max($x->[0], $x->[1]);

    my $start2 = min($y->[0], $y->[1]);
    my $end2   = max($y->[0], $y->[1]);

    if ($end1 >= $start2 && $end2 >= $start1){
        return min($end1, $end2) - max($start1, $start2)  + 1;
    }
    else {
        return 0;
    }
}


my %iupac = (
    G => [qw/ G /],
    A => [qw/ A /],
    T => [qw/ T /],
    C => [qw/ C /],
    R => [qw/ G A /],
    Y => [qw/ T C /],
    M => [qw/ A C /],
    K => [qw/ G T /],
    S => [qw/ G C /],
    W => [qw/ A T /],
    H => [qw/ A C T Y M W /],
    B => [qw/ G T C Y K S /],
    V => [qw/ G C A R M S /],
    D => [qw/ G A T R K W /],
    N => [qw/ G A T C R Y M K S W H B V D N/],
);

# base_match - return true if base1 is compatible with base2. 

sub base_match{
    my ($base1, $base2, $both) = @_;

    $base1 = uc $base1;
    $base2 = uc $base2;

    return first{ $base1 eq $_ } (@{$iupac{$base2}});
}

sub timestamp{ return strftime("%Y%m%d-%H%M",localtime); }
sub datestamp{ return strftime("%Y%m%d",localtime); }

sub fastq_read_length{
    my $filename = shift;
    open my $fh, "<", $filename;
    <$fh>;
    my $line = <$fh>;
    close $fh;

    if (defined $line){
        chomp $line;
        return length $line;
    }
    return;
}

sub read_conf{
    my ($index) = grep { '--conf' eq $main::ARGV[$_] } (0 .. $#main::ARGV);
    return () if ! defined $index;

    my (undef, $file) = splice @main::ARGV, $index, 2;
    return () if ! $file;

    die "--conf $file not a file?" if ! -f $file;

    return ParseConfig("$file");
}

1;


