#!/usr/bin/env perl

package My::SequenceWindower;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use DZUtil qw/safemethyl/;
use List::Util qw/min/;


sub new {
    my ($class, %opt) = @_;
    #my ($class, $sequence, $sequence_length, $feature, $window_size, $filename) = @_;
    my $self = bless {}, $class;

    $self->{sequence}        = $opt{sequence};
    $self->{sequence_length} = $opt{sequence_length};
    $self->{feature}         = $opt{feature};
    $self->{position}        = 1; # position should be aligned to windows
    $self->{last_position}   = 1; # last position to be seen via add_gff, for checking order
    $self->{window_size}     = $opt{window_size};
    $self->{no_skip}         = $opt{no_skip};

    $self->{c} = 0; # accumulators
    $self->{t} = 0;
    $self->{n} = 0;

    my $filename = $opt{filename};

    # file 
    if ($filename && -f $filename){
        open my $fh, '>', $filename;
        $self->{fh} = $fh;
    }
    else {
        $self->{fh} = \*STDOUT;
    }

    return $self;
}

sub in_current_window{
    my ($self, $new_position) = @_;
    my $pos = $self->{position}; 
    return $pos <= $new_position && $new_position < $pos + $self->{window_size};
}

sub add_gff{
    my ($self, $position, $c, $t, $n) = @_;
    if ($position < $self->{last_position}){
        die "file not ordered??";
    }
    elsif($self->in_current_window($position)){
        $self->{c} += $c;
        $self->{t} += $t;
        $self->{n} += $n;
    }
    # gff beyond current window
    else{
        # update position to next window
        $self->{position} += $self->{window_size};

        # fill_until the last window
        $self->fill_until($position);

        # initialize new window
        $self->{c} = $c;
        $self->{t} = $t;
        $self->{n} = $n;
    }
    $self->{last_position} = $position;
}

# fill from current position up to but not including window containing $until
sub fill_until{
    my ($self, $until) = @_;

    # with ws=50, round 299->250, 300->250, 301->300
    $until = ($until - 1) - ($until - 1) % $self->{window_size};

    while ($self->{position} <= $until){
        # if first window had mappings, blit that, and zero
        my ($c, $t, $n) = ($self->{c}, $self->{t}, $self->{n});
        if (grep { $_ > 0 } ($c, $t, $n)){
            # blit current
            say {$self->{fh}} join "\t", 
            $self->{sequence}, '.', $self->{feature}, 
            $self->{position}, $self->{position} + $self->{window_size} - 1,
            safemethyl($c, $t), ('.') x 2, "c=$c;t=$t;n=$n";
            ($self->{c}, $self->{t}, $self->{n}) = (0) x 3;
        }
        elsif ($self->{no_skip}){
            say {$self->{fh}}
            join "\t", 
            $self->{sequence}, 
            '.', 
            $self->{feature}, 
            $self->{position},
            min($self->{position} + $self->{window_size} - 1, $self->{sequence_length}),
            ('.') x 4;
        }
        $self->{position} += $self->{window_size};
    }
}

sub fill_tail{
    my $self = shift;
    $self->fill_until($self->{sequence_length} + $self->{window_size});
    # (add window size so any jagged end is included.
}

package main;
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;

END {close STDOUT}
$| = 1;

my $sw = My::SequenceWindower->new(
    sequence => 'chr1', 
    sequence_length => 59, 
    feature => 'CHH', 
    window_size => 20, 
    filename => undef,
    no_skip => 0,
);

#$sw->fill_until(301);
$sw->add_gff(3, 10, 20, 2);
$sw->add_gff(6, 10, 20, 2);
$sw->add_gff(7, 1, 7, 3);
$sw->add_gff(32, 1, 7, 3);
#, $position, $c, $t, $n) = @_;
$sw->fill_tail();


=head1 NAME

window_by_fixed_2.pl - even better window_by_fixed.pl
 
=head1 SYNOPSIS

 window_by_fixed_2.pl 

=head1 OPTIONS

=over

=item  -w <size> | --window-size <size>

=for Euclid
    size.default:     1
    size.type:        int, size >= 1 
    size.type.error:  <size> must be integer greater than 1

=item  -r <fasta> | --reference <fasta>

Reference fasta file.  If --no-skip/-k is used, for empty windows between the last window and the end of the
sequences will also be outputted.

=for Euclid
    fasta.type:        readable

=item -o <file> | --output <file>

=for Euclid 
    file.default: '-'

=item <files>...

=item -h | --help

=item -v | --verbose

=item  -n | --report-count 

Report 'n' in the scores column instead of c/(c+t).

=item --count-in-scores

Use the column 6 value instead of "n=#" for count score

=item  --debug <sqlite>

=for Euclid
    sqlite.type:        readable

=item  --keep-intermediate

=back

=cut

