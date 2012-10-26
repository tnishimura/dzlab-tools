package FastqReader::IsBS;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use FastqReader;
use DZUtil qw/c2t g2a reverse_complement open_filename_or_handle/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(is_bs);

sub is_bs{
    my $file = shift;
    my $fqr = FastqReader->new(file => $file);

    my ($a, $c, $g, $t) = (0) x 4;
    my $length;
    my $count = 1;
    while (defined(my $r = $fqr->next())){
        my ($readid, $sequence) = @$r;
        last if 10_000 == $count++;
        $length //= length $sequence;
        if (length $sequence != $length){
            die "$file: corrupt?";
        }
        $sequence = lc $sequence;

        $a += $sequence =~ tr/a/a/;
        $c += $sequence =~ tr/c/c/;
        $g += $sequence =~ tr/g/g/;
        $t += $sequence =~ tr/t/t/;
    }

    my $total = $a + $c + $g + $t;
    my $c_ratio = $c / $total;
    my $is_bs = $c_ratio < .125;

    if (wantarray){
        return ($is_bs, 
            total => $total, 
            a => $a/$total,
            c => $c/$total,
            g => $g/$total,
            t => $t/$total,
        );
    }
    else{
        return $is_bs;
    }

    #printf "%s: %.4f %s\n", $file, $c_ratio, ($c_ratio < .125 ? 'BS' : 'gDNA');
    #printf "a: %.4f\n", $a / $total;
    #printf "c: %.4f\n", $c / $total;
    #printf "g: %.4f\n", $g / $total;
    #printf "t: %.4f\n", $t / $total;
}

1;
