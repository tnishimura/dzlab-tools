package FastaReader::MethylSimulator;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Moose;
use Carp;
use autodie;    
use List::Util qw/sum/;

extends 'FastaReader';

# $ms->bsrecord()->{seq}{coord} = methylation count. 
# positive if positive strand, negative if reverse
has 'bsrecord' => (
    is => 'ro',
    default => sub { {} },
    init_arg => undef,
);

has 'read_length' => (
    is => 'ro',
    default => 100,
);

has 'methrate' => (
    is => 'ro',
    required => 1,
);

has 'norc' => (
    is => 'ro',
    default => 0, 
);

sub get_read{
    my $self = shift;

    my $read_length = $self->read_length();

    my %seqlen = $self->sequence_lengths();
    my $multargs = [(),()];

    while (my ($seq,$len) = each %seqlen) {
        push @{$multargs->[0]}, $seq;
        push @{$multargs->[1]}, $len;
    }

    my $seq = rmultinomial(@$multargs);

    my $start;
    my $stop;

    if ($read_length > $seqlen{$seq}){
        $start = 1;
        $stop = $seqlen{$seq};
    }
    else{
        $start = int(rand($seqlen{$seq}-$read_length))+1;
        $stop  = $start + $read_length -1;
    }

    my $get_rc   = $self->norc ? 0 : rand() > .50;

    my $read = $self->get($seq, $start, $stop, base => 1, rc => $get_rc);

    if ($self->methrate() > 0){
        my ($bsread, $meth) = bisulfite($read, $start, $stop, $get_rc, $self->methrate());

        $self->record_meth($seq, $get_rc, @$meth);

        return [$seq, $start, $stop, $get_rc, $read, $bsread];
    }
    else{
        return [$seq, $start, $stop, $get_rc, $read];
    }
}

sub record_meth{
    my ($self, $seq, $rc, @meth) = @_;
    my $bsrecord = $self->bsrecord();

    for my $pos (@meth) {
        # if forward strand, positive. if rc strand, negative. ok b/c there's no
        # way there's a C on both sides. I'm so crafty.
        if (($rc == 0 && exists $bsrecord->{$seq}{$pos} && $bsrecord->{$seq}{$pos} < 0) || 
            ($rc == 1 && exists $bsrecord->{$seq}{$pos} && $bsrecord->{$seq}{$pos} > 0)){
            croak "methylation on BOTH sides? BUG. $seq, $pos";
        }
        $bsrecord->{$seq}{$pos} += $rc == 0 ? 1 : -1;
    }
}

# $rate is methylation (protection from conversion)
# returns bisulfite simulated sequence and the absolute positions of methylation
sub bisulfite{
    my ($subseq, $start, $end, $rc, $rate) = @_;
    my @split = split //, $subseq;
    my @meth;
    for my $abspos ($start .. $end){
        my $relpos = $abspos - $start;
        if ($rc){
            $abspos = $end - $relpos;
        }

        if ($split[$relpos] eq 'C'){
            if( rand() > $rate){
                $split[$relpos] = 'T';
            }
            else{
                push @meth, $abspos;
            }
        }
    }
    @meth = sort { $a <=> $b } @meth;
    return join("", @split), \@meth; 
}

# rmultinomial [qw/A C G T/], [.25, .25, .25, .25]
sub rmultinomial{
    die "rmultinomial requires at least one probability-item pairs" if @_ < 2;
    my ($items, $probabilities) = @_;
    if (ref $items ne 'ARRAY' || ref $probabilities ne 'ARRAY' || @$items != @$probabilities){
        die "usage: rmultinomial [qw/A C G T/], [.25, .25, .25, .25]";
    }

    my $count = @$probabilities;
    my $norm= sum @$probabilities;
    @$probabilities = map { $_ / $norm } @$probabilities;

    if ($count == 1){
        return $items->[0];
    }

    my $unif = rand(1);
    
    my $cummulative_probability = 0;

    for my $i (0 .. $count - 1){

        $cummulative_probability += $probabilities->[$i];
        if ($unif < $cummulative_probability){
            return $items->[$i];
        }
    }
    
    # there's a chance we end up here b/c of floating point issues
    return $items->[-1];
}

sub dump{
    my ($self, $outfile) = @_;
    my %bsrec = %{$self->bsrecord()};
    open my $out, '>', $outfile;
    for my $seq (sort keys %bsrec) {
        for my $pos ( sort {$a <=> $b} keys %{$bsrec{$seq}} ) {
            my $score = $bsrec{$seq}{$pos};
            my $strand = $score > 0 ? '+' : '-';

            say $out join "\t", $seq, qw/. ./, $pos, $pos, abs($score), $strand, qw/. ./;
        }
    }
    close $out;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
