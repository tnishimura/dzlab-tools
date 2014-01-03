#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use List::MoreUtils qw/all any/;
use Term::ProgressBar;

use FindBin;
use lib "$FindBin::Bin/../lib";
use DZUtil qw/deprecation_message/;

deprecation_message("fastquack check");

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "fix|f=s" => \(my $fix),
    "max-error|e=i" => \(my $max_error = 1),
    "debug|d" => \(my $debug),
);
if (!$result){
    say "$0 [-f fixed.fastq] [-e 1] input.fastq";
    exit 1;
}

my $fix_fh;
if ($fix){
    open $fix_fh, '>', $fix;
}
my $nuc = qr/[ABCDGHKMNRSTVWY]+/i;
my $error_count = 0;

sub increment_error{
    my $error_msg = shift;
    if ($fix){
        if ($error_count < $max_error){
            $error_count++;
            say STDERR "WARNING: $error_msg at line $.: @_" if $debug;
        }
        else{
            die "error count exceeded, $error_msg at line $.: @_";
        }
    }
    else{
        die "$error_msg at line $.: @_";
    }
}

for my $file (@ARGV) {

    my $size = (stat($file))[7];

    say "$file:";

    my $pb = Term::ProgressBar->new({count => $size});
    my $counter = 0;
    $pb->minor(0);

    open my $in, '<', $file;

    LOOP:
    while (! eof $in){
        my $first_line = <$in>;
        if ($first_line =~ /\r\n$/){
            die "file has dos line ending, no soup for you.";
        }
        # check alignment of quartets.
        if ($first_line !~ /\@/){
            increment_error("quartet doesn't start with @", $first_line);
        }
        my @lines = ($first_line, map { scalar <$in> } (2 .. 4));
        my @lens = map { length $_ } @lines;

        # all undef
        if (all { ! defined $_ } @lines){
            last;
        }
        else{
            # only some undef
            if (grep { ! defined $_ } @lines){
                die "uneven number of lines @ around $."
            }
            elsif ($lines[0] !~ /^@/){
                increment_error("first line in quartet doesn't start with \@", @lines);
            }
            elsif ($lines[2] !~ /^\+/){
                increment_error("third line in quartet doesn't start with +", @lines);
            }
            elsif ($lens[0] != $lens[2] && $lines[2] !~ /^\+$/){
                increment_error("third line should either be a lone + or same length as line 1", @lines);
            }
            elsif ($lens[1] != $lens[3] ){
                increment_error("second and fourth line should be same length", @lines);
            }
            elsif ($lines[1] !~ m/^$nuc$/){
                increment_error("second should begin with $nuc", @lines);
            }
            elsif (any { m{[^[:ascii:]]} } @lines){
                increment_error("nonascii characters found", @lines);
            }
            elsif ($fix){
                print $fix_fh @lines;
            }
        }
        $pb->update(tell($in)) if ++$counter % 100_000 == 0;
    }
    $pb->update($size);

    close $in;
}
if ($fix) {
    close $fix_fh;
    say "OK, $error_count entries skipped";
}
else{
    say "OK";
}
exit 0;
