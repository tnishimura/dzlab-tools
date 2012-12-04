package GBUtil;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Bio::Graphics::Wiggle::Loader;
use File::Basename qw/basename/;
use Params::Validate qw/:all/;
use Hash::Util qw/lock_keys/;
use File::Which;
use Parallel::ForkManager;
use File::Path qw/make_path/;
use File::Spec::Functions qw/catfile rel2abs/;
use List::MoreUtils qw/notall/;
use YAML qw/LoadFile/;

use FastaReader;
use GFF::Parser;
use GFF::Statistics qw/gff_detect_width/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw( load_mysql_config prepare_gff prepare_fasta prepare_gff_to_wig );

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
            meta       => { default => undef, optional => 1, },
        });
    lock_keys(%opt);

    my $input_file_name = $opt{file};

    # setup file names
    croak "no such file $input_file_name" unless -f $input_file_name;
    my $staging_file_name = catfile($opt{stagingdir}, basename($input_file_name) . ".normalized.fasta");
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

# normalize gff. lowercase seqid. fill in source field (default to ".". this matches
# prepare_fasta's metafile's source).
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
    my $staging_file_name = catfile($opt{stagingdir}, basename($input_file_name) . ".normalized.gff");
    my %features;

    {
        open my $outfh, '>', $staging_file_name;
        my $p = GFF::Parser->new(file => $input_file_name);
        while (defined(my $gff = $p->next)){
            $gff->sequence(lc($gff->sequence() // '.'));
            $gff->source($opt{source} // '.');

            if (! exists $features{$gff->feature()}){
                $features{$gff->feature()} = 1;
            }

            say $outfh $gff;
        }
        close $outfh;
    }
    
    return ($staging_file_name, sort keys %features);
}

#######################################################################
# methylation gff

# compile wiggle to base_directory and write write meta file
sub compile_wiggle{
    my %opt = validate(@_, {
            feature        => 1, 
            source         => 1,
            track          => 1, # not sure what this is
            base_directory => 1, # where binary wiggles are placed
            meta_file      => 1,
            files => {
                type => ARRAYREF,
            },
        });
    lock_keys(%opt);

    my $loader = Bio::Graphics::Wiggle::Loader->new($opt{base_directory})
        or die "could not create loader";

    $loader->{trackname} = $opt{track} if defined $opt{track};

    for my $file (@{$opt{files}}) {
        say $file;
        my $fh = IO::File->new($file) or die "could not open $file: $!";
        $loader->load($fh);
    }
    open my $metafh, '>', $opt{meta_file};
    print {$metafh} $loader->featurefile('gff3',$opt{feature},$opt{source});
    close $metafh;
}

my @types = qw{methyl coverage};

sub _make_tmpname{
    my $file = shift;
    my $abs = rel2abs($file);
    return join "_", grep { $_ ne '' } split /\//, $abs;
}

#$meta_files{$feature}{$type} = sprintf($meta_filename_template, $feature, $type);
sub make_meta_filename{
    my %opt = validate(@_, { stagingdir => 1, file => 1, feature => 1, type => 1, });
    lock_keys(%opt);
    return sprintf(
        catfile($opt{stagingdir}, _make_tmpname($opt{file})) . "-%s-%s.meta.gff",
        $opt{feature}, $opt{type}
    );
}

# my $filename = sprintf($wig_filename_template, $seq, $feature, $type);
sub make_wig_filename{
    my %opt = validate(@_, { stagingdir => 1, file => 1, feature => 1, type => 1, seq => 1, });
    return sprintf(
        catfile($opt{stagingdir}, _make_tmpname($opt{file})) . "-%s-%s.%s.wig",
        $opt{seq}, $opt{feature}, $opt{type}
    );
}

# took out bed support b/c bed can support multiple seqids per file... should
# be in another subroutine.
# returns
#   [ 
#     { feature => 'CG', type => 'methyl | coverage', meta => 'meta.gff },
#     ...
#   ]
sub prepare_gff_to_wig{
    my %opt = validate(@_, {
            file       => 1,
            stagingdir => 1,
            wigdir     => 0,
            ctscore    => 0,
            source     => { default => undef, optional => 1, },
            compile    => { default => 1, optional     => 1, },
            parallel   => { default => 0, optional     => 1, },
        });
    lock_keys(%opt);

    make_path $opt{stagingdir};

    # name for staging files
    # my $wig_filename_template = catfile($opt{stagingdir}, basename($opt{file}, '.gff')) . "-%s-%s.%s.wig";
    # my $meta_filename_template = catfile($opt{stagingdir}, basename($opt{file}, '.gff')) . "-%s-%s.meta.gff";

    my $detected_width = gff_detect_width $opt{file};

    my %done_sequences;
    my %wig_files;  # $wig_files{$seq}{$feature}{methyl | coverage}{wig | handle} 
    my %meta_files; # $meta_files{$feature}{methyl | coverage}

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
        for my $type (@types) {
            if (!exists $meta_files{$feature}{$type}){
                $meta_files{$feature}{$type} = make_meta_filename(
                    file       => $opt{file},
                    stagingdir => $opt{stagingdir},
                    feature    => $feature,
                    type       => $type
                );
            }
            if (!exists $wig_files{$feature}{$type}{$seq}){
                my $filename = make_wig_filename(
                    file       => $opt{file},
                    stagingdir => $opt{stagingdir},
                    seq        => $seq,
                    feature    => $feature,
                    type       => $type
                );
                open my $fh, '>', $filename;
                $wig_files{$feature}{$type}{$seq}{name} = $filename;
                $wig_files{$feature}{$type}{$seq}{handle} = $fh;

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

        say {$wig_files{$feature}{methyl}{$seq}{handle}}   sprintf("%d\t%4f", $start, $score);
        say {$wig_files{$feature}{coverage}{$seq}{handle}} "$start\t", $c + $t;
    }

    # close handles and compile with compile_wiggle 
    for my $feature (keys %wig_files){
        for my $type (qw/methyl coverage/) {
            my @seqs = keys %{$wig_files{$feature}{$type}};

            close $wig_files{$feature}{$type}{$_}{handle} for @seqs;

            say STDERR "compiling $feature $type";

            compile_wiggle(
                files          => [map { $wig_files{$feature}{$type}{$_}{name} } @seqs],
                source         => $opt{source},
                track          => "$opt{source}-$feature-$type",
                feature        => "$feature-$type",
                base_directory => $opt{wigdir},
                meta_file      => $meta_files{$feature}{$type},
            );
        }
    }

    return map {
        my $feature = $_;
        map 
        { 
            {
                feature => $feature,
                type    => $_,
                meta    => $meta_files{$feature}{$_},
                source  => $opt{source},
            },
        } @types;
    } keys %meta_files;
}


1;
