package FastqReader::CheckIntegrity;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use List::MoreUtils qw/all any/;
use Term::ProgressBar;
use Params::Validate qw/:all/;
use IO::File;
use IO::Uncompress::Gunzip;
use IO::Uncompress::Bunzip2;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(check_integrity);

sub check_integrity{
    my %opt = validate(@_, {
            fix => {
                optional => 1,
            }, 
            maxerror => {
                type => SCALAR,
                regex => qr/^\d+$/,
                optional => 1,
            }, 
            files => {
                type => ARRAYREF,
                optional => 0,
            }, 
        });
    my $fix = $opt{fix};
    my $max_error = $opt{maxerror};
    my $debug = 1;
    my $files = $opt{files};

    my $nuc = qr/[ABCDGHKMNRSTVWY]+/i;

    FASTQFILE:
    for my $file (@$files) {
        if (!-f $file){
            croak "no such file $file";
        }
        say "**************************************************************";
        say "*** $file";
        say "";
        my $error_count = 0;
        my $increment_error = sub {
            my $error_msg = shift;
            if ($fix){
                if ($error_count < $max_error){
                    $error_count++;
                    say STDERR "WARNING: $error_msg at line $.:\n @_" if $debug;
                }
                else{
                    carp "$file: error count exceeded, $error_msg at line $.:\n @_";
                    next FASTQFILE;
                }
            }
            else{
                carp "$file: $error_msg at line $.: @_";
                next FASTQFILE;
            }
        };

        my $fix_fh;
        my $fix_file = "$file.FIXED";
        if ($fix){
            open $fix_fh, '>', $fix_file;
        }

        my $size = (stat($file))[7];

        my $pb = $file =~ /\.gz|bzip2$/ ? undef : Term::ProgressBar->new({count => $size});
        my $counter = 0;
        $pb->minor(0) if $pb;

        my $in = $file =~ /\.gz$/  ? IO::Uncompress::Gunzip->new($file) : 
                 $file =~ /\.bz2$/ ? IO::Uncompress::Bunzip2->new($file) : 
                 IO::File->new($file);

        LOOP:
        while (! eof $in){
            my $first_line = <$in>;
            if ($first_line =~ /\r\n$/){
                carp "$file: file has dos line ending, no soup for you.";
                next FASTQFILE;
            }
            # check alignment of quartets.
            if ($first_line !~ /\@/){
                $increment_error->("quartet doesn't start with @, skipping a line to see if simple extra/missing line", $first_line);
                next LOOP;
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
                    carp "$file: uneven number of lines @ around $.";
                    next FASTQFILE;
                }
                elsif ($lines[0] !~ /^@/){
                    $increment_error->("first line in quartet doesn't start with \@", @lines);
                }
                elsif ($lines[2] !~ /^\+/){
                    $increment_error->("third line in quartet doesn't start with +", @lines);
                }
                elsif ($lens[0] != $lens[2] && $lines[2] !~ /^\+$/){
                    $increment_error->("third line should either be a lone + or same length as line 1", @lines);
                }
                elsif ($lens[1] != $lens[3] ){
                    $increment_error->("second and fourth line should be same length", @lines);
                }
                elsif ($lines[1] !~ m/^$nuc$/){
                    $increment_error->("second should begin with $nuc", @lines);
                }
                elsif (any { m{[^[:ascii:]]} } @lines){
                    $increment_error->("nonascii characters found", @lines);
                }
                elsif ($fix){
                    print $fix_fh @lines;
                }
            }
            $pb->update(tell($in)) if ++$counter % 100_000 == 0 && $pb;
        }
        $pb->update($size) if $pb;

        close $in;
        if ($fix) {
            close $fix_fh;
            if ($error_count == 0){
                say "$file OK, no errors, fixing not necessary";
                unlink $fix_file;
            }
            else{
                say "$file had $error_count errors, fixed and created $fix_file";
            }
        }
        else{
            say "$file OK";
        }
    }
}

1;
