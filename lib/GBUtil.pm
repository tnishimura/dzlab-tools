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

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(gff_to_wig);

sub gff_to_wig{
    my %opt = validate(@_, {
            file    => 1,
            dir     => 0,
            compile => 0,
            bed     => 0,
            ctscore => 0,
        });

    my $file    = $opt{file};
    my $dir     = $opt{dir};
    my $compile = $opt{compile};
    my $bed     = $opt{bed};
    my $ctscore = $opt{ctscore};
    my $suffix  = $bed ? 'bed' : 'wig';

    # basename
    my $basename_template;
    if ($dir){
        $basename_template = catfile($dir, basename($file, '.gff')) . "-%s.%s.$suffix";
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
        my ($seq, $start, $end, $score, $c, $t) 
        = (lc($gff->sequence()), $gff->start(), $gff->end(), $gff->score(), $gff->get_column('c') // 0, $gff->get_column('t') // 0,);

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
            if (!exists $files{$seq}{$type}){
                my $filename = sprintf($basename_template, $seq, $type);
                open my $fh, '>', $filename;
                $files{$seq}{$type}{name} = $filename;
                $files{$seq}{$type}{handle} = $fh;

                say STDERR "created $filename";

                if (! $bed and defined $detected_width){
                    say $fh "variableStep  chrom=$seq  span=$detected_width";
                }
                elsif (! $bed and ! defined $detected_width){
                    say $fh "variableStep  chrom=$seq";
                }
                else{
                    say $fh qq{track type=wiggle_0 name="$filename $seq $type" };
                }
            }
        }

        if ($bed){
            say {$files{$seq}{methyl}{handle}}   sprintf("%d\t%d\t%4f", $start, $end, $score);
            say {$files{$seq}{coverage}{handle}} "$start\t$end\t", $c + $t;
        }
        else{
            say {$files{$seq}{methyl}{handle}}   sprintf("%d\t%4f", $start, $score);
            say {$files{$seq}{coverage}{handle}} "$start\t", $c + $t;
        }
    }

    # close and compile if necessary
    for my $seq (keys %files){
        for my $type (qw/methyl coverage/) {
            close $files{$seq}{$type}{handle};

            if ($compile){
                my $filename = $files{$seq}{$type}{name};
                say STDERR "compiling $filename";
                system("$wiggle2gff3 $filename > $filename.meta ");
            }
        }
    }
}


1;
