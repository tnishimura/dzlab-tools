#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Statistics qw/gff_info/;
use Template;
use Digest::MD5 qw/md5/;
use Pod::Usage;
use Getopt::Long;
use File::Spec::Functions qw/rel2abs/;


END {close STDOUT}
$| = 1;

#######################################################################
# stupid subroutine to assign colors 'randomly' but consistently,
# with defaults for common features.

{
    # colors that look ok on white.
    my @colors = qw{
        aqua 	aquamarine 	bisque 	black 	blue 	blueviolet 	brown 	burlywood 	cadetblue 	
        chartreuse 	chocolate 	coral 	cornflowerblue 	crimson 	cyan 	darkblue 	darkcyan 	
        darkgoldenrod 	darkgray 	darkgreen 	darkkhaki 	darkmagenta 	darkolivegreen 	darkorange 	darkorchid 	
        darkred 	darksalmon 	darkseagreen 	darkslateblue 	darkslategray 	darkturquoise 	darkviolet 	deeppink 	
        deepskyblue 	dimgray 	dodgerblue 	firebrick 	forestgreen 	fuchsia 	gold 	goldenrod 	
        gray 	green 	green 	greenyellow 	hotpink 	indianred 	indigo 	lawngreen 	
        lightblue 	lightcoral 	lightgreen 	lightpink 	lightsalmon 	lightseagreen 	lightskyblue 	lightslategray 	
        lightsteelblue 	lime 	limegreen 	magenta 	maroon 	mediumaquamarine 	mediumblue 	mediumorchid 	
        mediumpurple 	mediumseagreen 	mediumslateblue 	mediumslateblue 	mediumspringgreen 	mediumturquoise 	mediumvioletred 	midnightblue 	
        navy 	olive 	olivedrab 	orange 	orangered 	orchid 	palegreen 	paleturquoise 	
        palevioletred 	peachpuff 	peru 	pink 	plum 	powderblue 	purple 	red 	
        rosybrown 	royalblue 	saddlebrown 	salmon 	sandybrown 	seagreen 	sienna sienna 	
        silver skyblue slateblue slategray springgreen steelblue tan teal
        thistle tomato turquoise violet wheat yellow yellowgreen
    };

    # my %default_colors = (
    #     gene                 => "blue",
    #     exon                 => "purple",
    #     CDS                  => "red",
    #     chromosome           => "green",
    #     five_prime_UTR       => "brown",
    #     three_prime_UTR      => "brown",
    #     protein              => "cyan",
    #     mRNA                 => "yellow",
    #     transposable_element => "peachpuff",
    # );

    sub string2color{
        my $str = shift or die "need string";
        #if (exists $default_colors{$str}){ return $default_colors{$str}; }
        my ($index) = map { $_ % scalar(@colors) } unpack "l",  md5 $str;
        return $colors[$index];
    }
}

#######################################################################

my $result = GetOptions (
    "t|title=s"      => \my $title,
    "s|sqlite=s"     => \my $sqlite,
    "a|annotation=s" => \my $annotation,
    "g|gff=s"        => \my $gfffile,
    "o|output=s"     => \my $output,
);

if (!$result || ! $title || ! $sqlite){
    say "usage: gb-create-config.pl ...";
    exit 1;
}

#######################################################################

my @sequences;
my @features;
if ($annotation){
    my %annotation_info = %{gff_info $annotation};
    @sequences = keys %{$annotation_info{sequences}};
    @features = keys %{$annotation_info{features}};
}

#######################################################################
# parse gff for methyl

my @methyls;
my @coverages;
if ($gfffile){
    my $p = GFF::Parser->new(file => $gfffile, normalize => -1);

    while (defined(my $gff = $p->next)){
        my $feature = $gff->feature();
        if ($feature =~ /-methyl$/){
            push @methyls, $feature;
        }
        elsif ($feature =~ /-coverage$/){
            push @coverages, $feature;
        }
    }
    @methyls = sort @methyls;
    @coverages = sort @coverages;
}

#######################################################################
# TT

my $tt = Template->new() || die "$Template::ERROR\n";

my $vars = {
    title        => $title,
    sqlite       => rel2abs($sqlite),
    max_coverage => 1000,

    sequences    => \@sequences,
    features     => \@features,
    methyls      => \@methyls,
    coverages    => \@coverages,

    string2color => \&string2color,
};

# say Dumper $vars;

$tt->process(\*DATA, $vars, $output) || die $tt->error(), "\n";

__DATA__

[GENERAL]
description   = [% title %]
db_adaptor    = Bio::DB::SeqFeature::Store
db_args       = -adaptor DBI::SQLite
                -dsn     [% sqlite %]
plugins       = TrackDumper
examples      = [% FOREACH s = sequences %][% s %] [% END %]
                
# default features = Genes
# region segment = 10000
# initial landmark = chr1:5000..10000

########################
# Default glyph settings
########################

[TRACK DEFAULTS]
glyph         = generic
height        = 10
bgcolor       = lightgrey
fgcolor       = black
font2color    = blue
label density = 25
bump density  = 100
# where to link to when user clicks in detailed view
link          = AUTO

[DNA]
glyph          = dna
global feature = 1
height         = 40
do_gc          = 1
gc_window      = auto
fgcolor        = red
axis_color     = blue
strand         = both
category       = DNA
key            = DNA/GC Content

#######################################################################
# features found in annotation file

[% FOREACH f = features -%]
[[% f %]]
feature     = [% f %]
glyph       = generic
stranded    = 1
bgcolor     = [% string2color(f) %]
height      = 10
key         = [% f %]
description = 1

[% END %]

#######################################################################
# -methyl$ from gff files

[% FOREACH m = methyls -%]
[[% m %]]
feature        = [% m %]
glyph          = wiggle_xyplot
graph_type     = line
fgcolor        = black
bgcolor        = [% string2color(m) %]
height         = 50
min_score      = 0
max_score      = 1
scale          = right
group_on       = display_name
category       = Quantitative Data
key            = [% m %]

[% END %]

#######################################################################
# -coverage$ from gff files

[% FOREACH c = coverages -%]
[[% c %]]
feature        = [% c %]
glyph          = wiggle_xyplot
graph_type     = line
fgcolor        = black
bgcolor        = [% string2color(c) %]
height         = 50
min_score      = 0
max_score      = [% max_coverage %]
scale          = right
group_on       = display_name
category       = Quantitative Data
key            = [% c %]

[% END %]

