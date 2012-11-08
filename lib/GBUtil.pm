package GBUtil;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use File::Basename qw/basename/;
use Params::Validate qw/:all/;
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
our @EXPORT = qw(gff_to_wig load_mysql_config prepare_fasta prepare_gff);

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

# Bed output is probably not correct (at least gbrowse won't show it properly).
# leaving it in for now, but bed == GFF without binary compilation, so probably useless
# without fixing.
# return hash ($files{$seq}{$feature}{methyl | coverage}) 
sub prepare_gff_to_wig{
    my %opt = validate(@_, {
            file    => 1,
            dir     => 0,
            ctscore => 0,
            source  => 0, 
            compile => {
                default => 1,
                optional => 1,
            },
            bed => {
                default => 0,
                optional => 1,
            },
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
        $basename_template = ($file =~ s/(.*)\.gff$/$1/r) . "-%s.%s.%s.$suffix";
    }

    # wiggle2gff3 binary
    my $wiggle2gff3 = which('wiggle2gff3.pl') // which('wiggle2gff3') // undef;

    my $detected_width = gff_detect_width $file;

    my %done;
    my %files; # $files{$seq}{$feature}{methyl | coverage}{name | handle} 

    my $p = GFF::Parser->new(file => $file, normalize => 0);

    while (defined(my $gff = $p->next())){
        my ($seq, $start, $feature, $end, $score, $c, $t) 
        = (lc($gff->sequence()), $gff->start(), $gff->feature(), $gff->end(), $gff->score(), $gff->get_column('c') // 0, $gff->get_column('t') // 0,);

        my $width = $end - $start + 1;

        # okay for a single window to be different size b/c the last one may be smaller
        # if (defined $detected_width and $width != $detected_width){
        #     if (exists $done{$seq}){ die "uneven widths. line $.\n$gff"; }
        #     else{ $done{$seq} = 1; }
        # }

        $score //= "0";
        if ($ctscore and defined($c) and defined($t)){ $score = ($c + $t) == 0 ? 0 : $c / ($c + $t); }

        # create wigs if ! exist and write header
        for my $type (qw/methyl coverage/) {
            if (!exists $files{$seq}{$feature}{$type}){
                my $filename = sprintf($basename_template, $seq, $feature, $type);
                open my $fh, '>', $filename;
                $files{$seq}{$feature}{$type}{name} = $filename;
                $files{$seq}{$feature}{$type}{handle} = $fh;

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
            say {$files{$seq}{$feature}{methyl}{handle}}   sprintf("%s\t%d\t%d\t%4f", $seq, $start, $end, $score);
            say {$files{$seq}{$feature}{coverage}{handle}} sprintf("%s\t%d\t%d\t%d", $seq, $start, $end, $c + $t);
        }
        else{
            say {$files{$seq}{$feature}{methyl}{handle}}   sprintf("%d\t%4f", $start, $score);
            say {$files{$seq}{$feature}{coverage}{handle}} "$start\t", $c + $t;
        }
    }

    my $pm = Parallel::ForkManager->new($parallel);

    # close and compile if necessary
    for my $seq (keys %files){
        for my $feature (keys %{$files{$seq}}){
            for my $type (qw/methyl coverage/) {
                close $files{$seq}{$feature}{$type}{handle};
                delete $files{$seq}{$feature}{$type}{handle};

                if ($compile){
                    $pm->start and next;
                    my $filename = $files{$seq}{$feature}{$type}{name};
                    say STDERR "compiling $filename";
                    # system("$wiggle2gff3 --trackname=$feature-$type --method=$feature-$type $filename > $filename.meta ");
                    system("$wiggle2gff3 --method=$feature-$type $filename > $filename.meta ");
                    $pm->finish; 
                }
            }
        }
    }

    # rearrange hash to $files{$seq}{$feature}{methyl | coverage}
    for my $seq (keys %files){
        for my $feature (keys %{$files{$seq}}){
            for my $type (qw/methyl coverage/) {
                $files{$seq}{$feature}{$type} = $files{$seq}{$feature}{$type}{name};
            }
        }
    }
    
    $pm->wait_all_children;
    return %files;
}


# normalize fasta header lines, create associated GFF.
# my ($output_file_name, $meta_file_name) = prepare_fasta($input_file_name);
# perl -I$HOME/dzlab-tools/lib/ -MGBUtil -wle 'prepare_fasta("TAIR_reference.fa")'
sub prepare_fasta{
    my %opt = validate(@_, {
            file => {
                required => 1,
                callbacks => { 'exists' => sub { -f shift }, },
            },
            out => 0,
            meta => 0,
        });

    my $input_file_name = $opt{file};
    my $output_dir_or_file = $opt{out};

    croak "no such file $input_file_name" unless -f $input_file_name;
    my $output_file_name = 
       defined $output_dir_or_file && -d $output_dir_or_file ? 
          catfile($output_dir_or_file, basename($input_file_name) . ".normalized")
          : $input_file_name . ".normalized";
    my $meta_file_name = $opt{meta} // $output_file_name . ".meta";

    {
        open my $outfh, '>', $output_file_name;
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

# normalize gff. lowercase seqid. fill in source field.
# my $output_file_name = prepare_gff($input_file_name);
# perl -I$HOME/dzlab-tools/lib/ -MGBUtil -wle 'prepare_gff("foo.bar.gff")'
sub prepare_gff{
    my %opt = validate(@_, {
            file => {
                required => 1,
                callbacks => { 'exists' => sub { -f shift }, },
            },
            out => 0,
            source => 0,
        });

    my $input_file_name = $opt{file};
    my $output_dir_or_file = $opt{out};
    my $source = $opt{source} // basename($input_file_name);

    my $output_file_name = 
       defined $output_dir_or_file && -d $output_dir_or_file ? 
          catfile($output_dir_or_file, basename($input_file_name) . ".normalized")
          : $input_file_name . ".normalized";

    {
        open my $outfh, '>', $output_file_name;
        my $p = GFF::Parser->new(file => $input_file_name);
        while (defined(my $gff = $p->next)){
            $gff->sequence(lc($gff->sequence() // '.'));
            $gff->source($source);
            say $outfh $gff;
        }
        close $outfh;
    }
    
    return $output_file_name;
}

1;

__END__

parallel: 8
wiggledir: /tmp
fasta:
  - filename:
gff:
  - source: 
    filename: 
gffwig:
  - source: 
    filename: 
    ctscore:
