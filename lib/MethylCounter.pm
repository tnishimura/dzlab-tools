package MethylCounter;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;

use FastaReader;
use List::Util qw/max/;
use List::MoreUtils qw/all/;
use File::Temp qw/mktemp tempfile/;
use BigArray;

sub new {
    my $class = shift;
    my %opt = @_;
    my $self = bless {}, $class;

    $self->{bigarrays}        = {};
    $self->{stats}            = {};
    $self->{reference_genome} = $opt{reference_genome};
    $self->{dinucleotide}     = $opt{dinucleotide};
    $self->{file_prefix}      = $opt{file_prefix};

    return $self;
}


sub record_methylation{
    my ($self, $seq, $coord, $base) = @_;

    my $reference_genome = $self->{reference_genome};
    my $bigarrays = $self->{bigarrays};

    $seq = uc $seq;
    $base = uc $base;

    #say Dumper "BEFORE", $bigarrays;
    if (! exists $bigarrays->{$seq}){
        my $len = $reference_genome->get_length($seq);
        $bigarrays->{$seq}{C} = BigArray->new(size => $len, base => 1);
        $bigarrays->{$seq}{T} = BigArray->new(size => $len, base => 1);
    }
    #say Dumper "AFTER", $bigarrays;
    #say STDERR "$seq $base $coord";

    $bigarrays->{$seq}{$base}->push_increment($coord);
}

sub commit{
    my ($self) = @_;

    my $bigarrays = $self->{bigarrays};

    for my $seq (sort keys %$bigarrays) {
        $bigarrays->{$seq}{C}->commit_increment();
        $bigarrays->{$seq}{T}->commit_increment();
    }
}

sub record_output{
    my ($self) = @_;

    my $reference_genome = $self->{reference_genome};
    my $bigarrays        = $self->{bigarrays};
    my $stats            = $self->{stats};
    my $dinucleotide     = $self->{dinucleotide};
    my $file_prefix      = $self->{file_prefix};

    my @seqs = sort keys %$bigarrays;

    my %filehandles;
    my %temp2real;

    for my $seq (@seqs) {
        $bigarrays->{$seq}{C}->commit_increment();
        $bigarrays->{$seq}{T}->commit_increment();
    }

    for my $context (qw/CG CHG CHH/) {
        my $file = "$file_prefix.single-c.$context.gff.merged";
        my $tmpfile = mktemp($file . ".tmp.XXXX");
        open my $writer, q{>}, $tmpfile;
        $filehandles{$context} = $writer;
        $temp2real{$tmpfile} = $file;
    }

    %$stats = map {
        $_, 
        $dinucleotide ?  {
            unfiltered => { map {$_ => 0} qw/C CG CC CA CT T TG TC TA TT/ },
            filtered   => { map {$_ => 0} qw/C CG CC CA CT T TG TC TA TT/ },
        } : {
            unfiltered => { map {$_ => 0} qw/C CG CHH CHG T TG THH THG/ },
            filtered   => { map {$_ => 0} qw/C CG CHH CHG T TG THH THG/ },
        }
    } @seqs;

    for my $seq (@seqs) {
        my $len = $reference_genome->get_length($seq);

        my $cbigarray = $bigarrays->{$seq}{C};
        my $tbigarray = $bigarrays->{$seq}{T};

        POSITION:
        for my $coord (1 .. $len) {
            my $c = $cbigarray->get($coord);
            my $t = $tbigarray->get($coord);

            my $context = $reference_genome->get_context($seq, $coord, dinuc => $dinucleotide, undef_on_nonc => 1);
            next POSITION if ! defined $context; # means was not a C/G position

            $context = uc $context;

            die "why is \$c + \$t == 0? bug, dying" if $c + $t == 0;
            my $filtered = $c + $t > 20 ? '*' : '';

            my $score = sprintf("%.4f", $c/($c+$t));
            say {$filehandles{$context}} join "\t", $seq, q{.}, $context, $coord, $coord, $score, q{.}, q{.}, "c=$c;t=$t$filtered";

            my $type = $filtered ? 'filtered' : 'unfiltered';

            $stats->{$seq}{$type}{C} += $c;
            $stats->{$seq}{$type}{T} += $t;
            $stats->{$seq}{$type}{$context} += $c; # CG, CHG, CHH
            substr($context, 0, 1) = 'T';
            $stats->{$seq}{$type}{$context} += $t; # TG, THG, THH
        }
    }

    close $_ for values %filehandles;

    while (my ($tmp,$real) = each %temp2real) {
        rename $tmp, $real;
    }
}

sub count_methylation{
    my ($self, $correlation_gff) = @_;

    my $reference_genome = $self->{reference_genome};

    my ($seq, $start, $end, $reverse, $read_seq, $target_seq) = @$correlation_gff;
    $seq = uc $seq;
    my $filtered = 0;

    $self->{stats}{$seq}{bp} += length($read_seq);

    my @target_bases = split //,$target_seq;
    my @read_bases = (q{.}, q{.}, split(//, $read_seq), q{.}, q{.});

    READ:
    for (my $strand_coord = $start; $strand_coord <= $end; ++$strand_coord){
        my $i = $strand_coord - $start + 2;
        my $abs_coord = $reverse ? $reference_genome->get_length($seq) - $strand_coord + 1 : $strand_coord;
        my $context;

        my $read_base = $read_bases[$i];
        if ($target_bases[$i] eq 'C' && ($read_base eq 'C' || $read_base eq 'T')){
            $self->record_methylation($seq,$abs_coord,$read_base);
        }
        next READ;
    }
}

sub rat{
    my ($x,$y) = @_;
    if ($x+$y != 0){
        return sprintf("%.6f", $x/($x+$y));
    }
    else {
        return "0.000000";
    }
}


# this thing looks ridiculous...
sub print_freq{
    my ($self, $error) = @_;

    my $reference_genome = $self->{reference_genome};
    my $stats            = $self->{stats};
    my $dinucleotide     = $self->{dinucleotide};
    my $file_prefix      = $self->{file_prefix};

    open my $out, '>', "$file_prefix.freq";

    my @output;

    if ($dinucleotide){
        push @output, [qw/seq bp overlaps 
        C CG CT CA CC T TG TT TA TC
        C_ratio CG_ratio CT_ratio CA_ratio CC_ratio 
        filtered_C filtered_CG filtered_CT filtered_CA filtered_CC filtered_T filtered_TG filtered_TT filtered_TA filtered_TC
        filtered_C_ratio filtered_CG_ratio filtered_CT_ratio filtered_CA_ratio filtered_CC_ratio 
        /];
    }
    else{
        push @output, [qw/seq bp overlaps 
        C CG CHG CHH T TG THG THH 
        C_ratio CG_ratio CHG_ratio CHH_ratio 
        filtered_C filtered_CG filtered_CHG filtered_CHH filtered_T filtered_TG filtered_THG filtered_THH 
        filtered_C_ratio filtered_CG_ratio filtered_CHG_ratio filtered_CHH_ratio 
        /];
    }

    for my $seq (sort keys %$stats){
        if ($dinucleotide){
            push @output, [
            $seq, $stats->{$seq}{bp}, 0,

            $stats->{$seq}{unfiltered}{C}, $stats->{$seq}{unfiltered}{CG}, $stats->{$seq}{unfiltered}{CT},$stats->{$seq}{unfiltered}{CA}, $stats->{$seq}{unfiltered}{CC}, 
            $stats->{$seq}{unfiltered}{T}, $stats->{$seq}{unfiltered}{TG}, $stats->{$seq}{unfiltered}{TT},$stats->{$seq}{unfiltered}{TA}, $stats->{$seq}{unfiltered}{TC}, 

            rat($stats->{$seq}{unfiltered}{C}  ,$stats->{$seq}{unfiltered}{T} ),
            rat($stats->{$seq}{unfiltered}{CG} ,$stats->{$seq}{unfiltered}{TG} ),
            rat($stats->{$seq}{unfiltered}{CT} ,$stats->{$seq}{unfiltered}{TT} ),
            rat($stats->{$seq}{unfiltered}{CA} ,$stats->{$seq}{unfiltered}{TA} ),
            rat($stats->{$seq}{unfiltered}{CC} ,$stats->{$seq}{unfiltered}{TC} ),

            $stats->{$seq}{filtered}{C}, $stats->{$seq}{filtered}{CG}, $stats->{$seq}{filtered}{CT},$stats->{$seq}{filtered}{CA}, $stats->{$seq}{filtered}{CC}, 
            $stats->{$seq}{filtered}{T}, $stats->{$seq}{filtered}{TG}, $stats->{$seq}{filtered}{TT},$stats->{$seq}{filtered}{TA}, $stats->{$seq}{filtered}{TC}, 

            rat($stats->{$seq}{filtered}{C}  ,$stats->{$seq}{filtered}{T} ),
            rat($stats->{$seq}{filtered}{CG} ,$stats->{$seq}{filtered}{TG} ),
            rat($stats->{$seq}{filtered}{CT} ,$stats->{$seq}{filtered}{TT} ),
            rat($stats->{$seq}{filtered}{CA} ,$stats->{$seq}{filtered}{TA} ),
            rat($stats->{$seq}{filtered}{CC} ,$stats->{$seq}{filtered}{TC} ),
            ]
        }
        else{
            push @output, [
            $seq, $stats->{$seq}{bp}, 0,

            $stats->{$seq}{unfiltered}{C}, $stats->{$seq}{unfiltered}{CG}, $stats->{$seq}{unfiltered}{CHG}, $stats->{$seq}{unfiltered}{CHH},
            $stats->{$seq}{unfiltered}{T}, $stats->{$seq}{unfiltered}{TG}, $stats->{$seq}{unfiltered}{THG}, $stats->{$seq}{unfiltered}{THH},

            rat($stats->{$seq}{unfiltered}{C}  ,$stats->{$seq}{unfiltered}{T} ),
            rat($stats->{$seq}{unfiltered}{CG} ,$stats->{$seq}{unfiltered}{TG} ),
            rat($stats->{$seq}{unfiltered}{CHG},$stats->{$seq}{unfiltered}{THG} ),
            rat($stats->{$seq}{unfiltered}{CHH},$stats->{$seq}{unfiltered}{THH} ),

            $stats->{$seq}{filtered}{C}, $stats->{$seq}{filtered}{CG}, $stats->{$seq}{filtered}{CHG}, $stats->{$seq}{filtered}{CHH},
            $stats->{$seq}{filtered}{T}, $stats->{$seq}{filtered}{TG}, $stats->{$seq}{filtered}{THG}, $stats->{$seq}{filtered}{THH},

            rat($stats->{$seq}{filtered}{C}  , $stats->{$seq}{filtered}{T} ),
            rat($stats->{$seq}{filtered}{CG} , $stats->{$seq}{filtered}{TG} ),
            rat($stats->{$seq}{filtered}{CHG}, $stats->{$seq}{filtered}{THG} ),
            rat($stats->{$seq}{filtered}{CHH}, $stats->{$seq}{filtered}{THH} ),
            ]
        }
    }

    #say scalar(@$_) for @output;
    my $numcols = $dinucleotide ? 33 : 27;

    #die "uneven number of lines in freq? dying" unless all { $numcols == scalar @$_ } @output;

    for my $i (0..$numcols-1) {
        my $line = join "\t", map { $_->[$i] } @output;
        say STDERR $line;
        say $out $line;
    }
    say STDERR "error\t$error";
    say $out "error\t$error";
    close $out;

    #say Dumper \%stats;
}


1;

