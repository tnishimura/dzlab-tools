package GBUtil;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use File::Basename qw/basename/;
use GFF::Statistics qw/gff_detect_width/;
use Params::Validate qw/:all/;
use File::Which;
use Parallel::ForkManager;
use File::Path qw/make_path/;
use File::Spec::Functions qw/catfile/;
use List::MoreUtils qw/notall/;
use FastaReader;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(gff_to_wig load_mysql_config prepare_fasta);

# Bed output is probably not correct (at least gbrowse won't show it properly).
# leaving it in for now, but bed == GFF without binary compilation, so probably useless
# without fixing.
sub gff_to_wig{
    my %opt = validate(@_, {
            file    => 1,
            dir     => 0,
            compile => 0,
            bed     => 0,
            ctscore => 0,
            parallel => {
                default => 0,
                optional => 1,
            },
        });

    my $file    = $opt{file};
    my $dir     = $opt{dir};
    my $compile = $opt{compile};
    my $bed     = $opt{bed};
    my $ctscore = $opt{ctscore};
    my $parallel = $opt{parallel};
    my $suffix  = $bed ? 'bed' : 'wig';

    # basename
    my $basename_template;
    if ($dir){
        $basename_template = catfile($dir, basename($file, '.gff')) . "-%s-%s.%s.$suffix";
        make_path $dir;
    } 
    else{
        $basename_template = ($file =~ s/(.*)\.gff$/$1/r) . "-%s.%s.$suffix";
    }

    # wiggle2gff3 binary
    my $wiggle2gff3 = which('wiggle2gff3.pl') // which('wiggle2gff3') // undef;

    my $detected_width = gff_detect_width $file;

    my %done;
    my %files;

    my $p = GFF::Parser->new(file => $file, normalize => 0);

    while (defined(my $gff = $p->next())){
        my ($seq, $start, $feature, $end, $score, $c, $t) 
        = (lc($gff->sequence()), $gff->start(), $gff->feature(), $gff->end(), $gff->score(), $gff->get_column('c') // 0, $gff->get_column('t') // 0,);

        my $width = $end - $start + 1;

        # okay for a single window to be different size b/c the last one may be smaller
        if (defined $detected_width and $width != $detected_width){
            if (exists $done{$seq}){ die "uneven widths"; }
            else{ $done{$seq} = 1; }
        }

        $score //= "0";
        if ($ctscore and defined($c) and defined($t)){ $score = ($c + $t) == 0 ? 0 : $c / ($c + $t); }

        # create wigs if ! exist and write header
        for my $type (qw/methyl coverage/) {
            if (!exists $files{$seq}{$type}{$feature}){
                my $filename = sprintf($basename_template, $seq, $feature, $type);
                open my $fh, '>', $filename;
                $files{$seq}{$type}{$feature}{name} = $filename;
                $files{$seq}{$type}{$feature}{handle} = $fh;
                $files{$seq}{$type}{$feature}{feature} = "$feature-$type";

                say STDERR "created $filename";

                if (! $bed and defined $detected_width){
                    # say $fh qq{variableStep name="$feature-$type" chrom=$seq  span=$detected_width};
                    say $fh qq{variableStep chrom=$seq  span=$detected_width};
                }
                elsif (! $bed and ! defined $detected_width){
                    say $fh "variableStep  chrom=$seq";
                }
                else{
                    say $fh qq{track type=wiggle_0 name="Bed Format" description="BED format" visibility=full color=200,100,0 altColor=0,100,200 priority=20};
                }
            }
        }

        if ($bed){
            say {$files{$seq}{methyl}{$feature}{handle}}   sprintf("%s\t%d\t%d\t%4f", $seq, $start, $end, $score);
            say {$files{$seq}{coverage}{$feature}{handle}} sprintf("%s\t%d\t%d\t%d", $seq, $start, $end, $c + $t);
        }
        else{
            say {$files{$seq}{methyl}{$feature}{handle}}   sprintf("%d\t%4f", $start, $score);
            say {$files{$seq}{coverage}{$feature}{handle}} "$start\t", $c + $t;
        }
    }

    my $pm = Parallel::ForkManager->new($parallel);

    # close and compile if necessary
    for my $seq (keys %files){
        for my $type (qw/methyl coverage/) {
            for my $feature (keys %{$files{$seq}{$type}}){
                close $files{$seq}{$type}{$feature}{handle};

                if ($compile){
                    $pm->start and next;
                    my $filename = $files{$seq}{$type}{$feature}{name};
                    say STDERR "compiling $filename";
                    # system("$wiggle2gff3 --trackname=$feature-$type --method=$feature-$type $filename > $filename.meta ");
                    system("$wiggle2gff3 --method=$feature-$type $filename > $filename.meta ");
                    $pm->finish; 
                }
            }
        }
    }
    
    $pm->wait_all_children;
}

use YAML qw/LoadFile/;
# my ($user, $pass, $database, $host) = load_mysql_config();
sub load_mysql_config{
    my $file = shift // "$ENV{HOME}/.bioperl";
    croak "no such file $file" unless -f $file;

    my $config = LoadFile($file);
    # if (notall { exists $config->{$_} } qw/user pass database host/){
    #     croak "$file doesn't contain all of user, pass, database, host";
    # }
    return @{$config}{qw/user pass database host/};
}

# normalize fasta header lines, create associated GFF.
# my ($output_file_name, $meta_file_name) = prepare_fasta($input_file_name);
# perl -I$HOME/dzlab-tools/lib/ -MGBUtil -wle 'prepare_fasta("TAIR_reference.fa")'
sub prepare_fasta{
    my ($input_file_name, $output_file_name, $meta_file_name) = @_;
    croak "no such file $input_file_name" unless -f $input_file_name;
    $output_file_name //= $input_file_name . ".normalized";
    $meta_file_name //= $output_file_name . ".meta";

    {
        open my $outfh, '>', $output_file_name;
        open my $infh, '<:crlf', $input_file_name;
        while (defined(my $line = <$infh>)){
            chomp $line;
            if ($line =~ /^>(\w+)/){
                say $outfh ">$1";
            }
            else{
                say $outfh $line;
            }
        }
        close $infh;
        close $outfh;
    }
    
    my $fr = FastaReader->new(file => $output_file_name, normalize => 0, slurp => 0);
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

    return ($output_file_name, $meta_file_name);
}

1;
