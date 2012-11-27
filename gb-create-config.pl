#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Storable qw/dclone/;
use YAML qw/Load Dump LoadFile DumpFile/;

use FindBin;
use lib "$FindBin::Bin/lib";
use GBUtil;
use Digest::MD5 qw/md5/;

my $config_file = shift // die "need config file";
my $config = LoadFile($config_file);

# --- sample 
# fasta:
#   - meta: /home/toshiro/demeter/staging/TAIR_reference.fas.normalized.meta
#     source: TAIR8
#     staging: /home/toshiro/demeter/staging/TAIR_reference.fas.normalized
# gff:
#   - feature:
#       - CDS
#       - chromosome
#       - exon
#       - five_prime_UTR
#       - gene
#       - mRNA
#       - mRNA_TE_gene
#     source: TAIR8
#     staging: /home/toshiro/demeter/staging/TAIR8_gmod.gff.normalized
# gffwig:
#  - feature: CHG
#    meta: /home/toshiro/demeter/staging/all.chg-col.w50-CHG-methyl.meta.gff
#    source: at-en-lerfie-x-col-wt-chg
#    type: methyl
#  - feature: CHG
#    meta: /home/toshiro/demeter/staging/all.chg-col.w50-CHG-coverage.meta.gff
#    source: at-en-lerfie-x-col-wt-chg
#    type: coverage
#  - feature: CG
#    meta: /home/toshiro/demeter/staging/all.cg-col.w50-CG-methyl.meta.gff
#    source: at-en-lerfie-x-col-wt-cg
#    type: methyl
#  - feature: CG
#    meta: /home/toshiro/demeter/staging/all.cg-col.w50-CG-coverage.meta.gff
#    source: at-en-lerfie-x-col-wt-cg
#    type: coverage
#  - feature: CHH
#    meta: /home/toshiro/demeter/staging/all.chh-col.w50-CHH-methyl.meta.gff
#    source: at-en-lerfie-x-col-wt-chh
#    type: methyl
#  - feature: CHH
#    meta: /home/toshiro/demeter/staging/all.chh-col.w50-CHH-coverage.meta.gff
#    source: at-en-lerfie-x-col-wt-chh
#    type: coverage

use Tie::IxHash;
tie my %sections, 'Tie::IxHash'; 

my $database_nick = "scaffolds";

# GENERAL
say write_gbini_section( 
        GENERAL => {
            description  => "Desc",
            database     => $database_nick,
            plugins      => 'FilterTest RestrictionAnnotator TrackDumper FastaDumper',
            autocomplete => 1
            # 'initial landmark'  =>  'chr1:3000..6000',
            # default tracks = gene:TAIR8
        },

    );

# scaffolds:database
{
    my ($user, $pass, $database, $host) = load_mysql_config();
    say write_gbini_section( 
        "$database_nick:database" => {
            'db_adaptor'     =>  'Bio::DB::SeqFeature::Store',
            'db_args'        =>  qq{
            -adaptor DBI::mysql 
            -dsn dbi:mysql:database=$database;host=$host 
            -user $user 
            -pass $pass},
            # 'search options'  =>  ['all', 'default +autocomplete'],
            'search options'  =>  'default +autocomplete',
        }
    );
}

say write_gbini_section(
    'TRACK DEFAULTS' => {
        glyph           => 'generic',
        height          => 20,
        'label density' => 25,
        'bump density'  => 100,
        'show summary'  => 99999,  # go into summary mode when zoomed out to 100k
        database        => $database_nick,
        #bgcolor         => 'cyan',
        #fgcolor         => 'black',
    }
);


# GFF
for my $gff (sort { $a->{source} cmp $b->{source}} @{$config->{gff}}) {
    my $source = $gff->{source};
    for my $feature (@{$gff->{feature}}) {
        my $header = "$feature-$source";

        say write_gbini_section(
            $header => {
                category => $source,
                feature  => "$feature:$source",
                glyph    => 'gene',
                key      => "$feature $source",
                fgcolor  => string2color($header),
                height   => 10,
            }
        );
    }
}

# wiggles
for (sort 
    { 
        $a->{source}  cmp $b->{source} ||
        $a->{feature} cmp $b->{feature} ||
        $a->{type}    cmp $b->{type} 
    } @{$config->{gffwig}}){
    my ($feature, $meta, $source, $type) = @{$_}{qw/feature meta source type/};

    my $header = "$feature-$type-$source";
    my $body = {
        category     => $source,
        feature      => "$feature-$type:$source",
        glyph        => 'wiggle_xyplot',
        key          => "$feature $type $source",
        bgcolor      => string2color($header),
        graph_type   => 'linepoints',
        point_symbol => 'point',
        height       => 30,
        $type eq 'methyl' ? (
            max_score  =>  '1.0',
            min_score  =>  '0.0',
        ) :
        (),
    };

    say write_gbini_section($header => $body);
}

#######################################################################
# gb config file is not ini strictly, it allows multiple keys and multilines 
#         section_name => {
#           key1 => value1,
#           key2 => value2,
#           key3 => [value3a, value3b],
#         },

sub write_gbini_section{
    my ($name, $body) = @_;
    my @accum;
    if ($name ne '_'){
        push @accum, "[$name]"
    }
    for my $key (sort keys %$body) {
        my $value = $body->{$key};

        if (ref $value eq 'ARRAY'){
            for my $v (@$value) {
                $v = _indent_multilines($v);
                push @accum, "$key = $v";
            }
        }
        else{
            $value = _indent_multilines($value);
            push @accum, "$key = $value";
        }
    }
    push @accum, "";
    return join "\n", @accum;
}

# delete blank lines, indent second and subsequent lines
sub _indent_multilines{
    my $line = shift;
    return 
    join "\n    ", 
    grep { ! /^\s*$/ } 
    split /\n\s*/, 
    $line;
}

#######################################################################
# color assignment

# stupid subroutine to assign colors 'randomly' but consistently,
# with defaults for common features.

sub string2color{
    # colors that look ok on white.
    my @colors = qw{
    aqua aquamarine bisque black blue blueviolet brown burlywood cadetblue
    chartreuse chocolate coral cornflowerblue crimson cyan darkblue darkcyan
    darkgoldenrod darkgray darkgreen darkkhaki darkmagenta darkolivegreen
    darkorange darkorchid darkred darksalmon darkseagreen darkslateblue
    darkslategray darkturquoise darkviolet deeppink deepskyblue dimgray dodgerblue
    firebrick forestgreen fuchsia gold goldenrod gray green green greenyellow
    hotpink indianred indigo lawngreen lightblue lightcoral lightgreen lightpink
    lightsalmon lightseagreen lightskyblue lightslategray lightsteelblue lime
    limegreen magenta maroon mediumaquamarine mediumblue mediumorchid mediumpurple
    mediumseagreen mediumslateblue mediumslateblue mediumspringgreen
    mediumturquoise mediumvioletred midnightblue navy olive olivedrab orange
    orangered orchid palegreen paleturquoise palevioletred peachpuff peru pink plum
    powderblue purple red rosybrown royalblue saddlebrown salmon sandybrown
    seagreen sienna sienna silver skyblue slateblue slategray springgreen steelblue
    tan teal thistle tomato turquoise violet wheat yellow yellowgreen
    };
    my %default_colors = (
        gene                 => "blue",
        exon                 => "purple",
        CDS                  => "red",
        chromosome           => "green",
        five_prime_UTR       => "brown",
        three_prime_UTR      => "brown",
        protein              => "cyan",
        mRNA                 => "yellow",
        transposable_element => "peachpuff",
    );

    my $str = shift or die "need string";
    if (exists $default_colors{$str}){ return $default_colors{$str}; }
    my ($index) = map { $_ % scalar(@colors) } unpack "l",  md5 $str;
    return $colors[$index];
}

