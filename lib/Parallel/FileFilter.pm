
package Parallel::FileFilter;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(make_iterators);

sub make_iterators{
    my ($file,$num) = @_;
    $num ||= 1; # so you can pass 0, like p::fm.
    
    my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,
        $atime,$mtime,$ctime,$blksize,$blocks)
    = stat($file);

    my $step_size = int($size/$num);

    open my $fh, '<', $file;

    my $current_start = 0;
    my @accum;

    for my $i (1 .. $num) {
        my $current_end;
        if ($current_start + $step_size >=$size){
            $current_end = $size-1;
        }
        else{
            seek $fh,$step_size, 1;
            my $line = <$fh>;
            $current_end = tell $fh;
        }

        say STDERR "Creating iter for $current_start => $current_end (" . ($current_end - $current_start + 1) . ")";

        open my $closure_fh, '<', $file;
        seek $closure_fh, $current_start, 0;
        push @accum, sub {
            my $pos = tell $closure_fh;
            if ($pos != -1 && tell $closure_fh < $current_end){
                return scalar <$closure_fh>;
            }
            else {
                if ($pos != -1){
                    close $closure_fh;
                }
                return;
            }
        };
        $current_start = $current_end;
        if ($current_end == $size-1){
            last;
        }
    }
    close $fh;

    return @accum;
}

1;

