package MethylCounter;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use YAML qw/DumpFile/;

use FastaReader;
use List::Util qw/max/;
use List::MoreUtils qw/all/;
use File::Temp qw/mktemp tempfile/;
use BigArray;
use DZUtil qw/approximate_line_count/;
use GFF::Parser::Correlated;

use Params::Validate qw/:all/;
# types: SCALAR ARRAYREF HASHREF CODEREF GLOB GLOBREF SCALARREF UNDEF OBJECT(blessed) BOOLEAN(UNDEF | SCALAR) HANDLE

sub new {
    my $class = shift;
    my %opt = @_;
    my $self = bless {}, $class;

    $self->{dinucleotide} = $opt{dinucleotide};
    $self->{verbose}      = $opt{verbose};
    $self->{file}         = $opt{correlation};
    $self->{bigarrays}    = {};
    $self->{stats}        = {};

    $self->{genome} = FastaReader->new(
        file => $opt{genome}, 
        slurp => 1,
    );

    # get_length() memoized. 
    $self->{length} = {
        map {
            uc($_) => $self->{genome}->get_length($_)
        } $self->{genome}->sequence_list()
    };

    # want genome file instead of pre-created FastaReader b/c want to make sure to slurp

    return $self;
}

#######################################################################
# input

sub record_single_methylation{
    my ($self, $seq, $coord, $base) = @_;

    my $genome    = $self->{genome};
    my $bigarrays = $self->{bigarrays};

    $seq = uc $seq;
    $base = uc $base;

    if (! exists $bigarrays->{$seq}){
        #my $len = $genome->get_length($seq);
        my $len = $self->{length}{uc $seq};
        $bigarrays->{$seq}{C} = BigArray->new(size => $len, base => 1);
        $bigarrays->{$seq}{T} = BigArray->new(size => $len, base => 1);
    }

    $bigarrays->{$seq}{$base}->push_increment($coord);
}

sub count_methylation{
    my ($self, $correlation_gff) = @_;

    my $genome = $self->{genome};

    my ($seq, $start, $end, $reverse, $read_seq, $target_seq) = @$correlation_gff;
    $seq = uc $seq;
    my $filtered = 0;

    $self->{stats}{$seq}{bp} += length($read_seq);

    #my @target_bases = split //,$target_seq;
    #my @read_bases = (q{.}, q{.}, split(//, $read_seq), q{.}, q{.});
    $read_seq = "..$read_seq..";

    READ:
    for (my $strand_coord = $start; $strand_coord <= $end; ++$strand_coord){
        my $i = $strand_coord - $start + 2;
        # Reverse strand coordinates are w.r.t. the 3' end!!!! This is about the only place that this happens.
        #my $abs_coord = $reverse ? $genome->get_length($seq) - $strand_coord + 1 : $strand_coord; 
        my $abs_coord = $reverse ? $self->{length}{uc $seq} - $strand_coord + 1 : $strand_coord; 
        my $context;

        #my $read_base = $read_bases[$i];
        my $read_base = substr $read_seq, $i, 1;
        #if ($target_bases[$i] eq 'C' && ($read_base eq 'C' || $read_base eq 'T')){
        if (substr($target_seq, $i, 1) eq 'C' && ($read_base eq 'C' || $read_base eq 'T')){
            $self->record_single_methylation($seq,$abs_coord,$read_base);
        }
    }
}

sub process{
    my ($self) = @_;

    my $counter_increment = 10000;
    my $line_count        = approximate_line_count($self->{file}, 10000);
    my $verbose           = $self->{verbose};

    my $parser = GFF::Parser::Correlated->new(
        file => $self->{file}, 
        normalize => 0,
    );

    while (defined(my $corr = $parser->next())){
        if ($verbose && $. % $counter_increment == 0){
            printf(STDERR "%d (%.4f)\n", $., $. / $line_count * 100);
        }

        $self->count_methylation($corr);
    }
    say STDERR "Done processing! creating single-c and freq file" if $verbose;

    $self->{no_match_count} = $parser->no_match_counter();
    $self->{error_count} = $parser->error_counter();
}

#######################################################################
# output gff

sub commit_bigarrays {
    my ($self) = @_;

    my $bigarrays = $self->{bigarrays};

    for my $seq (sort keys %$bigarrays) {
        $bigarrays->{$seq}{C}->commit_increment();
        $bigarrays->{$seq}{T}->commit_increment();
    }
}

sub output_single_c { # and also count stats
    my ($self, %output_files) = @_;

    my $genome       = $self->{genome};
    my $bigarrays    = $self->{bigarrays};
    my $stats        = $self->{stats};
    my $dinucleotide = $self->{dinucleotide};

    my @contexts = $dinucleotide ? qw/CG CC CA CT/ : qw/CG CHG CHH/;

    croak '$mc->output_single(CG => cgfile, CHG => chgfile, CHH => chhfile)' 
    unless 3 == grep { exists $output_files{$_} } @contexts;

    my @seqs = sort keys %$bigarrays;

    my %filehandles;
    my %temp2real;

    $self->commit_bigarrays();

    for my $context (@contexts) {
        my $file = $output_files{$context};
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
        #my $len = $genome->get_length($seq);
        my $len = $self->{length}{uc $seq};

        my $cbigarray = $bigarrays->{$seq}{C};
        my $tbigarray = $bigarrays->{$seq}{T};

        POSITION:
        for my $coord (1 .. $len) {
            my $c = $cbigarray->get($coord);
            my $t = $tbigarray->get($coord);

            next POSITION if $c + $t == 0;

            my $filtered = $c + $t > 20 ? '*' : '';

            my $context = $genome->get_context($seq, $coord, dinuc => $dinucleotide, undef_on_nonc => 1);
            next POSITION if ! defined $context; # means was not a C/G position
            $context = uc $context;

            my $strand = uc($genome->get($seq, $coord, $coord)) eq 'G' ? '-' : '+';

            my $score = sprintf("%.4f", $c/($c+$t));

            say {$filehandles{$context}} join("\t", 
                $seq, q{.}, $context, $coord, $coord, $score, $strand, q{.}, "c=$c;t=$t$filtered",
            );

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

#######################################################################
# print_freq - must be run after output_single_c since stats need to be computed

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
    my ($self, $output_file) = @_;

    my $error = $self->{no_match_count} + $self->{error_count};
    my $genome       = $self->{genome};
    my $stats        = $self->{stats};
    my $dinucleotide = $self->{dinucleotide};

    my %out;
    for my $seq (sort keys %$stats){
        if ($dinucleotide){
            $out{$seq} = [
            { seq               => $seq},
            { bp                => $stats->{$seq}{bp}},
            { overlaps          => 0},
            { C                 => $stats->{$seq}{unfiltered}{C}},
            { CG                => $stats->{$seq}{unfiltered}{CG}},
            { CT                => $stats->{$seq}{unfiltered}{CT}},
            { CA                => $stats->{$seq}{unfiltered}{CA}},
            { CC                => $stats->{$seq}{unfiltered}{CC}},
            { T                 => $stats->{$seq}{unfiltered}{T}},
            { TG                => $stats->{$seq}{unfiltered}{TG}},
            { TT                => $stats->{$seq}{unfiltered}{TT}},
            { TA                => $stats->{$seq}{unfiltered}{TA}},
            { TC                => $stats->{$seq}{unfiltered}{TC}},
            { C_ratio           => rat($stats->{$seq}{unfiltered}{C}  ,$stats->{$seq}{unfiltered}{T} )},
            { CG_ratio          => rat($stats->{$seq}{unfiltered}{CG} ,$stats->{$seq}{unfiltered}{TG} )},
            { CT_ratio          => rat($stats->{$seq}{unfiltered}{CT} ,$stats->{$seq}{unfiltered}{TT} )},
            { CA_ratio          => rat($stats->{$seq}{unfiltered}{CA} ,$stats->{$seq}{unfiltered}{TA} )},
            { CC_ratio          => rat($stats->{$seq}{unfiltered}{CC} ,$stats->{$seq}{unfiltered}{TC} )},
            { filtered_C        => $stats->{$seq}{filtered}{C}},
            { filtered_CG       => $stats->{$seq}{filtered}{CG}},
            { filtered_CT       => $stats->{$seq}{filtered}{CT}},
            { filtered_CA       => $stats->{$seq}{filtered}{CA}},
            { filtered_CC       => $stats->{$seq}{filtered}{CC}},
            { filtered_T        => $stats->{$seq}{filtered}{T}},
            { filtered_TG       => $stats->{$seq}{filtered}{TG}},
            { filtered_TT       => $stats->{$seq}{filtered}{TT}},
            { filtered_TA       => $stats->{$seq}{filtered}{TA}},
            { filtered_TC       => $stats->{$seq}{filtered}{TC}},
            { filtered_C_ratio  => rat($stats->{$seq}{filtered}{C}  ,$stats->{$seq}{filtered}{T} )},
            { filtered_CG_ratio => rat($stats->{$seq}{filtered}{CG} ,$stats->{$seq}{filtered}{TG} )},
            { filtered_CT_ratio => rat($stats->{$seq}{filtered}{CT} ,$stats->{$seq}{filtered}{TT} )},
            { filtered_CA_ratio => rat($stats->{$seq}{filtered}{CA} ,$stats->{$seq}{filtered}{TA} )},
            { filtered_CC_ratio => rat($stats->{$seq}{filtered}{CC} ,$stats->{$seq}{filtered}{TC} )},
            ],
        }
        else {
            $out{$seq} = [
                { seq                => $seq},
                { bp                 => $stats->{$seq}{bp}},
                { overlaps           => 0},
                { C                  => $stats->{$seq}{unfiltered}{C}},
                { CG                 => $stats->{$seq}{unfiltered}{CG}},
                { CHG                => $stats->{$seq}{unfiltered}{CHG}},
                { CHH                => $stats->{$seq}{unfiltered}{CHH}},
                { T                  => $stats->{$seq}{unfiltered}{T}},
                { TG                 => $stats->{$seq}{unfiltered}{TG}},
                { THG                => $stats->{$seq}{unfiltered}{THG}},
                { THH                => $stats->{$seq}{unfiltered}{THH}},
                { C_ratio            => rat($stats->{$seq}{unfiltered}{C}  ,$stats->{$seq}{unfiltered}{T} )},
                { CG_ratio           => rat($stats->{$seq}{unfiltered}{CG} ,$stats->{$seq}{unfiltered}{TG} )},
                { CHG_ratio          => rat($stats->{$seq}{unfiltered}{CHG},$stats->{$seq}{unfiltered}{THG} )},
                { CHH_ratio          => rat($stats->{$seq}{unfiltered}{CHH},$stats->{$seq}{unfiltered}{THH} )},
                { filtered_C         => $stats->{$seq}{filtered}{C}},
                { filtered_CG        => $stats->{$seq}{filtered}{CG}},
                { filtered_CHG       => $stats->{$seq}{filtered}{CHG}},
                { filtered_CHH       => $stats->{$seq}{filtered}{CHH}},
                { filtered_T         => $stats->{$seq}{filtered}{T}},
                { filtered_TG        => $stats->{$seq}{filtered}{TG}},
                { filtered_THG       => $stats->{$seq}{filtered}{THG}},
                { filtered_THH       => $stats->{$seq}{filtered}{THH}},
                { filtered_C_ratio   => rat($stats->{$seq}{filtered}{C}  , $stats->{$seq}{filtered}{T} )},
                { filtered_CG_ratio  => rat($stats->{$seq}{filtered}{CG} , $stats->{$seq}{filtered}{TG} )},
                { filtered_CHG_ratio => rat($stats->{$seq}{filtered}{CHG}, $stats->{$seq}{filtered}{THG} )},
                { filtered_CHH_ratio => rat($stats->{$seq}{filtered}{CHH}, $stats->{$seq}{filtered}{THH} )},
            ]
        }
    }

    DumpFile $output_file, \%out;
}

#######################################################################
# run

sub run{
    use Params::Validate qw/:all/;
    # types: SCALAR ARRAYREF HASHREF CODEREF GLOB GLOBREF SCALARREF UNDEF OBJECT(blessed) BOOLEAN(UNDEF | SCALAR) HANDLE
    
    my %opt = validate(@_, {
            dinucleotide => 0, 
            genome => {
                isa => 'FastaReader',
            },
            file_prefix => 1,
        });

}


1;

