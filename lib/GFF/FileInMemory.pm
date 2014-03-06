package GFF::FileInMemory;
use v5.12.0;
use Moose;
use Moose::Exporter;
use Data::Dumper;
use Carp;
use GFF::Parser;
use autodie;    

has file => (is => 'ro', required => 1);
has normalize => (is => 'ro', required => 0, default => 0);
has ignore_bad_coords => (is => 'ro', default => 1);
has gffs => (
    traits  => ['Array'],
    is      => 'ro',
    isa     => 'ArrayRef[GFF]',
    default => sub { [] },
    handles => {
        add_gff       => 'push',
        map_gffs      => 'map',
        filter_gffs   => 'grep',
        find_gff      => 'first',
        get_by_index => 'get',
        has_gffs      => 'count',
        count         => 'count',
        sort_in_place => 'sort_in_place',
    },
);

sub _sorter_triple{
    my ($seq_a, $seq_b, $start_a, $start_b, $end_a, $end_b) = @_;
    $seq_a   cmp $seq_b   || $start_a <=> $start_b || $end_a   <=> $end_b 
}

sub _sorter_double{
    my ($seq_a, $start_a, $seq_b, $start_b) = @_;
    return ($seq_a cmp $seq_b)   || ($start_a <=> $start_b );
}

sub sort_by_position{
    my $self = shift;
    $self->sort_in_place(sub { 
            _sorter_triple(
                $_[0]->sequence , $_[0]->sequence ,
                $_[0]->start    , $_[0]->start    ,
                $_[0]->end      , $_[0]->end
            )
        });
}

sub print_to_file{
    my ($self, $file) = @_;
    open my $fh, '>', $file;
    for my $gff (@{$self->gffs}) {
        say $fh $gff;
    }
    close $fh;
}

sub sequences_in_file{
    my ($self) = @_;
    my %seqs = $self->map_gffs(sub { defined($_->sequence) ? ($_->sequence => 1) : () });
    return sort keys %seqs;
}

sub find_by_seq_start{
    my ($self, $target_seq, $target_start) = @_;

    my $gffs = $self->gffs;

    my $min_index = 0;
    my $max_index = $self->count - 1;

    while ($max_index >= $min_index){
        my $mid_point = int(($max_index + $min_index) / 2);
        my $try = $gffs->[$mid_point];
        my $try_seq = $try->sequence;
        my $try_start = $try->start;
        my $cmp = _sorter_double($try->sequence, $try->start, $target_seq, $target_start,);
        # say "($min_index, $mid_point, $max_index) $cmp = _sorter_double($try_seq, $try_start, $target_seq, $target_start)";
        if ($cmp== 0){
            return $try;
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

sub BUILD{
    my ($self) = @_;
    my $parser = GFF::Parser->new(file => $self->file, normalize => $self->normalize);
    while (defined(my $gff = $parser->next)){
        if ($self->ignore_bad_coords && (! defined $gff->sequence || ! defined $gff->start)){
            next;
        }
        $self->add_gff($gff);
    }
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;

