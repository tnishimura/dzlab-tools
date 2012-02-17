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


require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(base_composition methylation_context);

sub make_iterator{
    my ($fr, $seq, $motif_length) = @_;
    my @buffer;
    my $position = 0;
    my $seqlen = $fr->get_length($seq);
    my $bufsize = 1000;

    return sub {
        if (@buffer < $motif_length){
            if ($position < $seqlen){
                push @buffer, split //, $fr->get($seq, $position, $position + $bufsize - 1, base => 0, lenient => 1);
                $position += $bufsize;
            }
            else{
                return;
            }
        }
        my @motif = @buffer[0 .. $motif_length - 1];
        shift @buffer;
        return \@motif
    };
}

sub base_composition{
    my ($file, $motif_length) = @_;
    my $fr = FastaReader->new(file => $file, slurp => 0);

    my %score;

    # prepopulate, so that 0 counts are not omitted
    my $base_iterator = _base_iterator(4, $motif_length);
    my @bases = qw/A C T G/;
    while (defined(my $n = $base_iterator->())){
        for my $seq ($fr->sequence_list()) {
            $score{$seq}{join "", map { $bases[$_] } @$n} = 0;
        }
    }

    for my $seq ($fr->sequence_list()) {
        my $iter = make_iterator($fr, $seq, $motif_length);
        while (defined(my $motif = $iter->())){
            my $joined = uc join "", @$motif;
            if ($joined =~ /^[ACGT]+$/){
                ++$score{$seq}{join "", @$motif};
            }
        }
    }
    return \%score;
}

use List::Util qw/sum/;

sub methylation_context{
    my ($file,$bc3) = @_;
    $bc3 //= base_composition $file, 3;

    my %scores;

    my $fr = FastaReader->new(file => $file, slurp => 0);
    for my $seq ($fr->sequence_list()) {
        $scores{$seq}{CHG} = sum map { $bc3->{$seq}{$_} } qw/CAG CCG CTG CGG/;
        $scores{$seq}{CHH} = sum map { $bc3->{$seq}{$_} } qw/CAA CCA CTA CGA CAC CCC CTC CGC CAT CCT CTT CGT/;
    }
    return \%scores;
}


sub report{
    my ($file, $max) = @_;

    my @bc;
    $bc[0] = base_composition $file, 1;
    $bc[1] = base_composition $file, 2;
    $bc[2] = base_composition $file, 3;
    for (4 .. $max){
        $bc[$_ - 1] = base_composition $file, $_;
    }

    my $meth = methylation_context($file, $bc[2]);

    say Dumper \@bc;
    my @seq_list = sort keys %{$bc[0]};

    say join "\t", "Context", @seq_list;

    for my $scores ($bc[0], $meth, @bc[1 .. $max]){
        for my $context (sort keys %{$scores->{$seq_list[0]}}){
            say join "\t", $context, map { $scores->{$_}{$context} } (@seq_list);
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

sub bc_table{
    

}

1;
