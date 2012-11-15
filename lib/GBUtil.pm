package GBUtil;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use File::Basename qw/basename/;
use Params::Validate qw/:all/;
use Hash::Util qw/lock_keys/;
use File::Which;
use Parallel::ForkManager;
use File::Path qw/make_path/;
use File::Spec::Functions qw/catfile/;
use List::MoreUtils qw/notall/;
use YAML qw/LoadFile/;

use FastaReader;
use GFF::Parser;
use GFF::Statistics qw/gff_detect_width/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw( load_mysql_config prepare_gff_to_wig prepare_gff prepare_fasta );

# my ($user, $pass, $database, $host) = load_mysql_config();
sub load_mysql_config{
    my $file = shift // "$ENV{HOME}/.bioperl";
    croak "no such file $file" unless -f $file;

    my $config = LoadFile($file);
    croak "can't read $file, malformatted?" unless ref $config eq 'HASH';
    # if (notall { exists $config->{$_} } qw/user pass database host/){
    #     croak "$file doesn't contain all of user, pass, database, host";
    # }
    return @{$config}{qw/user pass database host/};
}

# <old>Bed output is probably not correct (at least gbrowse won't show it properly).
# leaving it in for now, but bed == GFF without binary compilation, so probably useless
# without fixing.</old> 
# took out bed support b/c bed can support multiple seqids per file... should
# be in another subroutine.
# return hash ($meta_files{$seq}{$feature}{methyl | coverage}) 
sub prepare_gff_to_wig{
    my %opt = validate(@_, {
            file     => 1,
            stagingdir => 1,
            wigdir   => 0,
            ctscore  => 0,
            source   => { default => undef, optional => 1, },
            compile  => { default => 1, optional => 1, },
            parallel => { default => 0, optional => 1, },
        });
    lock_keys(%opt);

    # name for staging files
    my $staging_filename_template = catfile($opt{stagingdir}, basename($opt{file}, '.gff')) . "-%s-%s.%s.wig";
    make_path $opt{stagingdir};

    # wiggle2gff3 binary
    my $wiggle2gff3 = which('wiggle2gff3.pl') // which('wiggle2gff3') // undef;

    # binary wiggle output dir
    $opt{wigdir} //= $opt{stagingdir};
    my $wiggle_path = defined $opt{wigdir} ? "--path=$opt{wigdir}" : "";

    my $detected_width = gff_detect_width $opt{file};

    my %done_sequences;
    my %staging_files; # $staging_files{$seq}{$feature}{methyl | coverage}{name | handle} 
    my %meta_files;    # $meta_files{$seq}{$feature}{methyl | coverage}

    my $p = GFF::Parser->new(file => $opt{file}, normalize => 0);

    # convert gff to wig 
    while (defined(my $gff = $p->next())){
        my ($seq, $start, $feature, $end, $score, $c, $t) 
        = (lc($gff->sequence()), $gff->start(), $gff->feature(), $gff->end(), $gff->score(), $gff->get_column('c') // 0, $gff->get_column('t') // 0,);

        my $width = $end - $start + 1;

        # okay for a single window to be different size b/c the last one may be smaller
        if (defined $detected_width and $width != $detected_width){
            if (exists $done_sequences{$seq}){ die "uneven widths. line $.\n$gff"; }
            else{ $done_sequences{$seq} = 1; }
        }

        $score //= "0";
        if ($opt{ctscore} and defined($c) and defined($t)){ $score = ($c + $t) == 0 ? 0 : $c / ($c + $t); }

        # create wigs if ! exist and write header
        for my $type (qw/methyl coverage/) {
            if (!exists $staging_files{$seq}{$feature}{$type}){
                my $filename = sprintf($staging_filename_template, $seq, $feature, $type);
                open my $fh, '>', $filename;
                $staging_files{$seq}{$feature}{$type}{name} = $filename;
                $meta_files{$seq}{$feature}{$type} = "$filename.meta";
                $staging_files{$seq}{$feature}{$type}{handle} = $fh;

                say STDERR "created $filename";

                if (defined $detected_width){
                    # say $fh qq{variableStep name="$feature-$type" chrom=$seq  span=$detected_width};
                    say $fh qq{variableStep chrom=$seq  span=$detected_width};
                }
                elsif (! defined $detected_width){
                    say $fh "variableStep  chrom=$seq";
                }
            }
        }

        say {$staging_files{$seq}{$feature}{methyl}{handle}}   sprintf("%d\t%4f", $start, $score);
        say {$staging_files{$seq}{$feature}{coverage}{handle}} "$start\t", $c + $t;
    }

    my $pm = Parallel::ForkManager->new($opt{parallel});

    # close handles and compile with wiggle2gff3 if necessary
    for my $seq (keys %staging_files){
        for my $feature (keys %{$staging_files{$seq}}){
            for my $type (qw/methyl coverage/) {
                close $staging_files{$seq}{$feature}{$type}{handle};
                delete $staging_files{$seq}{$feature}{$type}{handle};

                if ($opt{compile}){
                    $pm->start and next;
                    my $staging = $staging_files{$seq}{$feature}{$type}{name};
                    my $meta = $meta_files{$seq}{$feature}{$type};
                    say STDERR "compiling $staging";
                    # system("$wiggle2gff3 --trackname=$feature-$type --method=$feature-$type $filename > $filename.meta ");
                    system("$wiggle2gff3 $wiggle_path --method=$feature-$type $staging > $meta");
                    $pm->finish; 
                }
            }
        }
    }

    $pm->wait_all_children;
    return %meta_files;
}

# normalize fasta header lines, create associated GFF.
# my ($staging_file_name, $meta_file_name) = prepare_fasta($input_file_name);
# perl -I$HOME/dzlab-tools/lib/ -MGBUtil -wle 'prepare_fasta("TAIR_reference.fa")'
sub prepare_fasta{
    my %opt = validate(@_, {
            file => {
                required => 1,
                callbacks => { 'exists' => sub { -f shift }, },
            },
            stagingdir => 1,
            meta     => { default => undef, optional => 1, },
        });
    lock_keys(%opt);

    my $input_file_name = $opt{file};

    # setup file names
    croak "no such file $input_file_name" unless -f $input_file_name;
    my $staging_file_name = catfile($opt{stagingdir}, basename($input_file_name) . ".normalized");
    my $meta_file_name = $opt{meta} // $staging_file_name . ".meta";

    say "+++ $staging_file_name";

    # copy input to meta, normalizing headers in the process.
    {
        open my $outfh, '>', $staging_file_name;
        open my $infh, '<:crlf', $input_file_name;
        while (defined(my $line = <$infh>)){
            chomp $line;
            if ($line =~ /^>(\w+)/){
                say $outfh ">\L$1";
            }
            else{
                say $outfh $line;
            }
        }
        close $infh;
        close $outfh;
    }
    
    # create meta
    my $fr = FastaReader->new(file => $staging_file_name, normalize => 0, slurp => 0);
    my %seqlen = $fr->sequence_lengths;

    open my $metafh, '>', $meta_file_name;
    say $metafh "##gff-version 3\n";
    for my $seqid (sort keys %seqlen) {
        say $metafh join "\t", 
        $seqid, 
        qw/./, # ? 
        'chromosome',
        1,
        $seqlen{$seqid},
        qw/. . ./,
        "ID=$seqid;Name=$seqid",
    }
    close $metafh;

    return ($staging_file_name, $meta_file_name);
}

# normalize gff. lowercase seqid. fill in source field (default to basename of input file).
# my $staging_file_name = prepare_gff($input_file_name);
# perl -I$HOME/dzlab-tools/lib/ -MGBUtil -wle 'prepare_gff("foo.bar.gff")'
sub prepare_gff{
    my %opt = validate(@_, {
            file => {
                required => 1,
                callbacks => { 'exists' => sub { -f shift }, },
            },
            stagingdir => 1,
            source => {
                default => undef,
                optional => 1,
            },
        });
    lock_keys(%opt);

    my $input_file_name = $opt{file};
    my $staging_file_name = catfile($opt{stagingdir}, basename($input_file_name) . ".normalized");

    my $source = $opt{source} // basename($input_file_name, '.gff');

    {
        open my $outfh, '>', $staging_file_name;
        my $p = GFF::Parser->new(file => $input_file_name);
        while (defined(my $gff = $p->next)){
            $gff->sequence(lc($gff->sequence() // '.'));
            $gff->source($source);
            say $outfh $gff;
        }
        close $outfh;
    }
    
    return $staging_file_name;
}

1;

__END__
# example

parallel:    4
stagingdir:  /home/toshiro/demeter/staging
wigdir:      /home/toshiro/demeter/wig
fasta:
  - file:    /home/toshiro/genomes/AT/TAIR_reference.fas
gff:
  - file:    /home/toshiro/annotations/AT/gmod/TAIR8_gmod.gff
gffwig:
  - file:    /home/toshiro/GEO-Submission-AT-Demeter-2012/windows/at-endosperm-ler_fie_x_col_wt/windows-Col/all.cg-col.w50.gff
    source:  at-en-lerfie-x-col-wt-cg
    ctscore: 0 
  - file:    /home/toshiro/GEO-Submission-AT-Demeter-2012/windows/at-endosperm-ler_fie_x_col_wt/windows-Col/all.chg-col.w50.gff
    source:  at-en-lerfie-x-col-wt-chg
    ctscore: 0 
  - file:    /home/toshiro/GEO-Submission-AT-Demeter-2012/windows/at-endosperm-ler_fie_x_col_wt/windows-Col/all.chh-col.w50.gff
    source:  at-en-lerfie-x-col-wt-chh
    ctscore: 0 
