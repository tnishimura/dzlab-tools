package Ends::NeighborMapCollection;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Moose;
use Carp;
use autodie;    
use Ends::NeighborMap;
use GFF::Parser;
use Tree::Range;
use overload '""' => \&stringify;

has file            => (is => 'ro', required => 1,);
has tag             => (is => 'ro', required => 1, );
has flag            => (is  => 'ro', required => 1, );
has distance        => (is  => 'ro', required => 1, );
has prime           => (is  => 'ro', required => 1, );
has flag_6_distance => (is  => 'ro', required => 1, );
has binwidth        => (is  => 'ro', required => 1, );
has numbins         => (is  => 'rw', );


# id => [start_bin, end_bin]
# keep track of which bins are possible.
# also, all_id grabs id list.
has valid_bin_range => (
    traits    => ['Hash'],
    is        => 'ro',
    isa       => 'HashRef[Any]',
    default   => sub { {} },
    handles   => {
        set_valid_bin_range     => 'set',
        get_valid_bin_range     => 'get',
        all_id   => 'keys',  # BUG: only keeps track of id's with bins.  overlapped genes won't appear, need to handle separately
    },
);

# hash of sequence => Tree::Range of [id, bin#]
has lookup_tree => (
    traits    => ['Hash'],
    is        => 'ro',
    default   => sub { {} },
    handles   => {
        set_lookup_tree => 'set',
        get_lookup_tree => 'get',
        has_lookup_tree => 'exists',
        all_seq   => 'keys'
    },
);

has nm => (
    traits    => ['Hash'],
    is        => 'ro',
    isa       => 'HashRef[Any]',
    default   => sub { {} },
    handles   => {
        set_nm     => 'set',
        get_nm     => 'get',
        list_nm   => 'keys'
    },
);
sub lookup{
    my ($self, $seq, $start, $end) = @_;
    my $tree = $self->get_lookup_tree(uc $seq);
    return $tree->search_overlap($start,$end);
}
sub bin_valid{
    my ($self, $id, $bin) = @_;
    my ($start, $end) = @{$self->get_valid_bin_range($id)};
    return ($start <= $bin && $bin <= $end);
}

sub stringify{
    my $self = shift;
    for my $kv (sort $self->list_nm) {
        say $self->get_nm($kv);
    }
}

sub BUILD{
    my ($self) = @_;

    my $binwidth = $self->binwidth;
    my $distance = $self->distance;
    my $flag = $self->flag;
    my $flag_6_distance = $self->flag_6_distance;
    my $prime = $self->prime;
    $self->numbins(my $numbins = int(($self->distance * 2) / $binwidth));

    my %neighbormaps;
    my %lookup_trees;
    my %id_list;

    my $parser = GFF::Parser->new(file => $self->file);

    # load gff file into neighbormaps
    while (defined(my $gff = $parser->next())){
        my ($seq, $start, $end, $strand, $id) = 
        ($gff->sequence, $gff->start, $gff->end, $gff->strand, $gff->get_column($self->tag));
        next if ! defined $seq; 
        $seq = uc $seq;
        push @{$id_list{$seq}}, $id;
        if (! exists $neighbormaps{$seq}){
            $neighbormaps{$seq} = 
            Ends::NeighborMap->new(
                flag            => $flag,
                distance        => $distance,
                prime           => $prime,
                flag_6_distance => $flag_6_distance,
            );
        }
        $neighbormaps{$seq}->add($id, $strand, $start, $end);
    }

    # finalize neighbor maps, load into trees
    for my $seq (keys %id_list) {
        my $nm = $neighbormaps{$seq};
        my $tr = Tree::Range->new();
        $nm->finalize();

        IDLOOP:
        for my $id (@{$id_list{$seq}}) {
            my ($position, $strand, $overlapped, $flag_upstream, $flag_downstream)
            = $nm->neighborhood($id);

            next IDLOOP if $overlapped;
            
            my $binnum = $strand eq '+' ? 0 : $numbins - 1; # 0->99 or 99->0
            my $first_valid_bin = $numbins;
            my $last_valid_bin = 0;
            my $start = $position - $distance + 1;
            my $end = $start + $binwidth - 1;
            # -299->-200, -199->-100, -99->0, 1->100, 101->200, 201->300
            
            while ($start <= $distance + $position){
                if (Tree::Range::_overlap($start, $end, $flag_upstream, $flag_downstream)){
                    $tr->add($start, $end, [$id, $binnum]); 
                    $first_valid_bin = $binnum < $first_valid_bin ? $binnum : $first_valid_bin;
                    $last_valid_bin = $last_valid_bin < $binnum ? $binnum : $last_valid_bin;
                }
                $binnum += $strand eq '+' ? 1 : -1;
                $start += $binwidth;
                $end += $binwidth;
            }
            if ($first_valid_bin < $last_valid_bin){
                $self->set_valid_bin_range($id, [$first_valid_bin, $last_valid_bin]);
            }
            else {
                $self->set_valid_bin_range($id, [-1, -1]);
            }
        }
        $tr->finalize();
        $self->set_lookup_tree($seq,$tr);
        $self->set_nm($seq,$nm);
    }
}

use Storable;
sub cache_name{
    my %args = @_;
    sprintf "%s.armacache-%s.%s.%s.%s.%s.%s",
    $args{file},           
    $args{tag},            
    $args{flag},           
    $args{distance},       
    $args{prime},          
    $args{flag_6_distance},
    $args{binwidth},       
}

sub new_cached{
    my %args = @_;
    my $cache = cache_name(%args);
    if (-f $cache){
        say STDERR "armageddon cache found at $cache, using that instead of rebuilding";
        my $ret = retrieve($cache);
        if (defined $ret){
            return $ret;
        }
        else{
            say STDERR "cache seems broken, will delete and rebuild";
            if (! unlink $cache){
                die "couldn't even delete... something wrong";
            }
        }
    }
    say STDERR "armageddon will cache at $cache, for later";
    my $nmc = Ends::NeighborMapCollection->new(%args);
    store $nmc, $cache;
    return $nmc;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
