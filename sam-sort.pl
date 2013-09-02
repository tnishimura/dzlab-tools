#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;

use File::Temp;
use File::Copy;

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "by-readid|r" => \(my $by_readid),
    "by-position|p" => \(my $by_position),
);
my $sam = shift;
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || ! ($by_readid xor $by_position) || ! $sam) ;  

my $unsorted = "$sam.unsorted";
my $tmp = "$sam.tmp";

if (-f $unsorted){
    die "$unsorted exists?";
}

sub cleanup{
    if (-f $unsorted){
        move($unsorted, $sam);
    }
    if (-f $tmp){
        unlink $tmp;
    }
    exit;
}

$SIG{INT} = sub{
    cleanup;
    exit 1;
};

END{
    cleanup;
}

move($sam, $unsorted);

# grab header lines, rest into tmp
my @header_lines;
{
    my $tmpfh = IO::File->new($tmp, "w");
    my $unsorted_fh = IO::File->new($unsorted);
    while (defined(my $line = <$unsorted_fh>)){
        if ($line =~ /^@/){
            push @header_lines, $line;
        }
        else{
            $tmpfh->print($line);
        }
    }
    $unsorted_fh->close;
    $tmpfh->close;
}

# sort

if ($by_readid){
    if (0 != system("sort", "-S10%", "-k1,1", "-o", $tmp, $tmp)){
        die "sort failed?";
    }
}
elsif ($by_position){
    if (0 != system("sort", "-S10%", "-k2,2", "-k3,3n", "-o", $tmp, $tmp)){
        die "sort failed?";
    }
}

# copy saved header and sorted tmp to original sam file
{
    my $tmpfh = IO::File->new($tmp);
    my $samfh = IO::File->new($sam, "w");

    $samfh->print($_) for @header_lines;

    while (defined(my $line = <$tmpfh>)){
        $samfh->print($line);
    }
    $samfh->close;
    $tmpfh->close;
    unlink $unsorted;
    unlink $tmp;
}

=head1 sam-sort.pl 

Usage examples:

 sam-sort.pl {--by-readid|--by-position} input.sam
 sam-sort.pl {-r|-p} input.sam

=cut

