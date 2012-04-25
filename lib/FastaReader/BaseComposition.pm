package FastaReader::BaseComposition;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use FastaReader;
use Tie::IxHash;
use List::Util qw/sum/;
use DZUtil qw/reverse_complement/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(base_composition);

our $VERBOSE = 0;

sub make_iterator{
    my ($fr, $seq, $motif_length) = @_;
    my @buffer;
    my $position = 0;
    my $seqlen = $fr->get_length($seq);
    my $bufsize = 50000;

    return sub {
        if (@buffer < $motif_length){
            if ($position < $seqlen){
                push @buffer, split //, $fr->get($seq, $position, $position + $bufsize - 1, base => 0, lenient => 1);
                $position += $bufsize;
            }
        }
        if (@buffer == 0){
            return;
        }
        my @motif = @buffer[($motif_length <= @buffer) ? (0  .. $motif_length - 1) : (0 .. $#buffer)];
        shift @buffer;
        return \@motif
    };
}

sub base_composition{
    my ($file, $motif_length, $methyl) = @_;
    my $size = [stat($file)]->[7];
    my $fr = FastaReader->new(file => $file, slurp => $size < 2**30);

    croak "if you want to do methylation scores, length needs to be at least 3"
    if $methyl && $motif_length < 3;

    my %composite_score;        # everything
    my %rc_score;
    my %single_score; # slightly redundant... 


    # prepopulate, so that 0 counts are not omitted. A C G T AA AC AG AT .. TG TT ..
    my @bases = qw/A C T G/;
    for (1 .. $motif_length){
        my $base_iterator = _base_iterator(4, $_);
        while (defined(my $n = $base_iterator->())){
            for my $seq ($fr->sequence_list()) {
                my $context = join "", map { $bases[$_] } @$n;

                if (length $context == 1){
                    $single_score{$seq}{$context} = 0;
                }
                else{
                    $composite_score{$seq}{$context} = 0;
                    # for CHG, CHH, we need to record both strands
                    if ($methyl && ($context eq 'CG' || $context =~ /^C..$/)){
                        $rc_score{$seq}{$context} = 0;
                    }
                }
            }
        }
    }
    #die Dumper \%score, \%rc_score;

    for my $seq ($fr->sequence_list()) {
        my $iter = make_iterator($fr, $seq, $motif_length);
        my $counter = 0;

        while (defined(my $motif = $iter->())){
            for (0 .. scalar(@$motif) -1){
                my $joined = uc join "", @{$motif}[0 .. $_];
                if ($joined =~ /^[ACGT]+$/){
                    if (length $joined == 1){
                        ++$single_score{$seq}{$joined};
                    }
                    else{
                        ++$composite_score{$seq}{$joined};
                        $joined = reverse_complement $joined;
                        if ($methyl && ($joined eq 'CG' || $joined =~ /^C..$/)){
                            ++$rc_score{$seq}{$joined};
                        }
                    }

                }
            }
            say STDERR "$seq $counter/" . $fr->get_length($seq) if ++$counter % 100_000 == 0 && $VERBOSE;
        }

        # methylation count if >=3 
        if ($methyl){
            $composite_score{$seq} = {
                CG => ($composite_score{$seq}{CG} + $rc_score{$seq}{CG}),
                CHG => (sum map { $composite_score{$seq}{$_} + $rc_score{$seq}{$_} } qw/CAG CCG CTG/),
                CHH => (sum map { $composite_score{$seq}{$_} + $rc_score{$seq}{$_} } qw/CAA CAC CAT CCA CCC CCT CTA CTC CTT/),
            };
        }
    }

    return (\%single_score, \%composite_score);
}


sub report{
    my ($file, $max, $methylation) = @_;

    my ($single, $compo) = base_composition $file, $max, $methylation;
    my @seq_list = sort keys %$single;

    #say Dumper $score;
    say join "\t", "Context", @seq_list;

    for my $context (qw/A C G T/){
        say join "\t", $context, map { $single->{$_}{$context} } @seq_list;
    }

    my @contexts = sort keys %{$compo->{$seq_list[0]}};
    for my $l (1 .. $max) {
        for my $context ( grep { length($_) == $l } @contexts){
            say join "\t", $context, map { $compo->{$_}{$context} } @seq_list;
        }
    }
}


# convert a base10 number to digits in base $base, wrapping around.
# convert_base(10,2,4) => [1,0,1,0]
# convert_base(65536,2,4) => [0,0,0,0]
sub _convert_base {
    my ($num, $base, $wraparound_digits) = @_;
    my @accum = ();

    while ($num > 0){
        my $remainder = $num % $base;
        push @accum, $remainder;
        $num = int(0.5 + ($num - $remainder)/$base);
    }
    if (@accum <= $wraparound_digits){
        return [reverse map { $accum[$_] // 0 } (0 .. $wraparound_digits - 1)];
    }
    else{
        return [reverse @accum[0..$wraparound_digits-1]];
    }
}

# return [0, 0, 0, 0], [0, 0, 0, 1], .., [1,1,1,1], ... wraparound
sub _base_iterator{
    my ($base, $wraparound_digits) = @_;
    my $counter = 0;
    sub {
        return if ($counter >= $base ** $wraparound_digits);
        return _convert_base($counter++, $base, $wraparound_digits);
    }
}

1;