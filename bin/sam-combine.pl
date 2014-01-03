#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;
use List::MoreUtils qw/uniq/;

use FindBin;
use lib "$FindBin::Bin/../lib";

use Sam::Parser;
use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "output|o=s" => \(my $output = '-'),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result);  

my $outfh = $output eq '-' ? \*STDOUT : IO::File->new($output, 'w');

# uniq from List::MoreUtils is stable
my @common_headers = uniq map 
{
    # use Sam::Parser only for getting headers. when getting body, 
    # just grep out ^@ since it's faster
    my $p = Sam::Parser->new(file => $_);
    @{$p->header_lines};
} @ARGV;

for my $h (@common_headers) {
    $outfh->print("$h\n");
}

for my $f (@ARGV) {
    my $in = IO::File->new($f);
    while (defined(my $line = <$in>)){
        $outfh->print($line) if $line !~ /^@/;
    }
}

$outfh->close;

=head1 sam-combine.pl 

Combine sam files, taking care of headers.

Usage examples:

 sam-combine.pl -o out.sam sam1.sam sam2.sam ...

=cut

