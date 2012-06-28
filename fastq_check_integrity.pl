#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use List::MoreUtils qw/all/;
use Term::ProgressBar;

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "fix|f=s" => \(my $fix),
    "max-error|e=i" => \(my $max_error = 1),
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
        # check alignment of quartets.
        if ($first_line !~ /\@/){
            if ($fix){
                if ($error_count < $max_error){
                    $error_count++;
                    next LOOP;
                }
                else{
                    die "error count exceeded, quartet doesn't start with @ at line: $.";
                }
            }
            else{
                die "quartet doesn't start with @ at line: $.";
            }
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
            elsif ($lines[0] !~ /^@/ 
                || $lines[2] !~ /^\+/ 
                # third line could be just a lone '+', in newer fastq files
                || ($lens[0] != $lens[2] && $lines[2] !~ /^\+$/)
                || $lens[1] != $lens[3] 
                || $lines[1] !~ m/^$nuc$/
            ){
                if ($fix){
                    if ($error_count < $max_error){
                        $error_count++;
                        next LOOP;
                    }
                    else{
                        die "error count exceeded, last quartet was @lines";
                    }
                }
                else{
                    die "malformed FASTQ quad @ $.\n" . join "", @lines;
                }
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
