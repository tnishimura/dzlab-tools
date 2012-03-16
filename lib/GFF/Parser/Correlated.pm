package GFF::Parser::Correlated;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Moose;
use Carp;
use autodie;    

extends 'GFF::Parser';

has 'correlation_error_counter' => (
    traits  => ['Counter'],
    is      => 'ro',
    isa     => 'Num',
    default => 0,
    handles => {
        inc_correlation_error_counter   => 'inc',
    },
);

has 'no_match_counter' => (
    traits  => ['Counter'],
    is      => 'ro',
    isa     => 'Num',
    default => 0,
    handles => {
        inc_no_match_counter   => 'inc',
    },
);

# die on # unrecoverable error
# return (0, $errmsg) on continuable error
#        (1, [$seq, $start, $end, $filetered, $reverse, $read_seq, $target_seq]) on 
around 'next' => sub{
    my $super = shift;
    my $self = shift;

    LOOP:
    while (defined(my $gff = $self->$super())){
        # body...
        my ($seq, undef, $col3, $start,$end,undef, $strand, undef, $col9) = @{$gff->as_array()};

        my $filtered = (defined $col9 && $col9 =~ s/\*$//) ? 1 : 0;

        if (!defined $seq){
            $self->inc_no_match_counter();
            next LOOP;
        }

        # parse for read, target sequence

        my ($read_seq, $target_seq);
        if ($col3 =~ /([A-Z]+)$/){
            $read_seq= $1;
        }
        if ($col9 =~ /target=([A-Z]+)$/){
            $target_seq = $1;
        }

        croak "strand should be + or -" unless $strand ~~ [qw/+ -/];
        croak "can't get read sequence from column 3" if ! defined $read_seq;
        croak "can't get target sequence from column 3" if ! defined $target_seq;
        croak "read seq is not the right length" if (length($read_seq) != ($end-$start+1));
        croak "target seq is not length(read)+4" if (length $target_seq != 4 + length $read_seq);

        if (! only_c2t_changes($read_seq, $target_seq)){
            $self->inc_correlation_error_counter();
            next LOOP;
        }

        my $reverse = $strand eq '+' ? 0 : $strand eq '-' ? 1 : croak "strand should be + or -";

        return [$seq, $start, $end, $filtered, $reverse, $read_seq, $target_seq];
    }
    return;

};

# VERY temporary kludgely hack to avoid errors in correlation.gff...
# returns false if more than a quarter are non c2t changes...
sub only_c2t_changes{
    my ($read_str, $target_str) = @_;
    $read_str =~ tr/C/T/;
    $target_str =~ tr/C/T/;
    my @read = split //, $read_str;
    my @target = split //, $target_str;

    my $total = 0;
    for my $index (0 .. $#read) {
        my ($r, $t) = ($read[$index], $target[2+$index]);
        if ($r ne $t){
            $total++;
        }
    }

    return ($total/length $read_str) < 0.25 ? 1 : 0;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
