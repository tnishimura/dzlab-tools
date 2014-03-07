package GFF::DataFrame;
use v5.12.0;
use Moose;
use Moose::Exporter;
use Data::Dumper;
use Carp;
use autodie;    

use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;

package GFF::DataFrame::Coord{
    use v5.12.0;
    use warnings FATAL => "all";
    use autodie;
    use base qw(Class::Accessor);
    __PACKAGE__->mk_accessors(qw(sequence start end));
    use overload q("") => sub { 
        my $self = shift;
        return join ";", $self->sequence, $self->start, $self->end // ".";
    };

    sub new_from_gff{
        my (undef, $gff) = @_;
        GFF::DataFrame::Coord->new({
                sequence => $gff->sequence,
                start    => $gff->start,
                end      => $gff->end,
            }
        );
    }

    # class methods
    sub compare{
        my (undef, $left, $right) = @_;
        return (lc($left->sequence) cmp lc($right->sequence)) || 
               ($left->start        <=> $right->start) || 
               ($left->end          <=> $right->end);
    }
    sub compare_no_end{
        my (undef, $left, $right) = @_;
        return (lc($left->sequence) cmp lc($right->sequence)) || 
               ($left->start        <=> $right->start);
    }
    sub find_coord_index{
        my (undef, $coords, $target_sequence, $target_start, $target_end) = @_;
        my $target_coord = GFF::DataFrame::Coord->new( {sequence => $target_sequence, start => $target_start, end => $target_end});

        my $min_index = 0;
        my $max_index = scalar(@$coords) - 1;

        while ($max_index >= $min_index){
            my $mid_point = int(($max_index + $min_index) / 2);
            my $try = $coords->[$mid_point];
            my $cmp = defined $target_end 
                        ? GFF::DataFrame::Coord->compare($try, $target_coord)
                        : GFF::DataFrame::Coord->compare_no_end($try, $target_coord);
                        #say "($min_index, $mid_point, $max_index) $cmp = compare($try, $target_coord)";
            if ($cmp== 0){
                return $mid_point;
            }
            elsif ($cmp < 0){ # try coord is less than target coord
                # say "setting min_index to " . ($mid_point + 1);
                $min_index = $mid_point + 1;
            }
            elsif ($cmp > 0){ # try coord is greater than target coord
                # say "setting max_index to " . ($mid_point - 1);
                $max_index = $mid_point - 1;
            }
            else{
                die "doubleyouteaeff?"
            }
        }
        return;
    }

}

sub new_from_gffs{
    shift;
    my @files = @_;
    croak "need files" if ! @files;
    my %coords_accum;
    my %cells_accum;
    # first pass - read in coords into hash in order to get uniques,
    # put values into cells_accum to transform later
    for my $f (@files) {
        my $p = GFF::Parser->new(file => $f);
        while (defined(my $gff = $p->next)){
            if (! defined $gff->sequence ||
                ! defined $gff->start    || 
                ! defined $gff->end){
                next;
            }
            my $coord = GFF::DataFrame::Coord->new_from_gff($gff);
            $coords_accum{"$coord"} = $coord;
            push @{$cells_accum{$f}}, $gff;
        }
    }

    my @coords = sort { GFF::DataFrame::Coord->compare($a, $b) } values %coords_accum;
    my $num_coords = @coords;
    my %cells;

    for my $f (@files) {
        $cells{$f} = [(undef) x $num_coords];
        for my $gff (@{$cells_accum{$f}}) {
            my $index = GFF::DataFrame::Coord->find_coord_index(\@coords, $gff->sequence, $gff->start, $gff->end);
            die "impossible for gff to not have matching coord here..." if (! defined $index);

            $cells{$f}[$index] = $gff;
        }
    }

    return __PACKAGE__->new(
        coords => \@coords,
        cells => \%cells,
    );
}

has coords => (is => 'ro', required => 1);
has cells  => (is => 'ro', required => 1);

sub files           { my $self = shift; return keys %{$self->cells}; }
sub num_coordinates { my $self = shift; return scalar(@{$self->coords}); }

sub get_gffs{
    my ($self, $seq, $start, $end) = @_;
    my $index = GFF::DataFrame::Coord->find_coord_index($self->coords, $seq, $start, $end);
    return map { $_ => $self->cells->{$_}[$index] } $self->files;
}

sub make_iterator{
    my $self = shift;
    my @indices = (0 .. $self->num_coordinates - 1);
    sub {
        my $index = shift @indices;
        if (defined $index){
            my $coord = $self->coords->[$index];
            my %files2gffs = map { $_ => $self->cells->{$_}[$index] } $self->files;
            return (
                $coord->sequence,
                $coord->start,
                $coord->end,
                %files2gffs,
            );
        }
        return ();
    };
}

sub dump{
    my $self = shift;
    my $iter = $self->make_iterator;
    while (my @results = $iter->()){
        my ($sequence, $start, $end, %files2gffs) = @results;
        say "$sequence, $start, $end";
        for my $f (sort keys %files2gffs) {
            say $files2gffs{$f} // ".";
        }
    }
}

# my $x = __PACKAGE__->new_from_gffs( "t/data/table1.gff", "t/data/table2.gff",);
# $x->dump;

no Moose;
__PACKAGE__->meta->make_immutable;

1;
