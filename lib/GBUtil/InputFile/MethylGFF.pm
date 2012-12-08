package GBUtil::InputFile::MethylGFF;
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
use File::Path qw/make_path/;
use File::Spec::Functions qw/catfile rel2abs/;
use Readonly;
use GFF::Parser;
use GFF::Statistics qw/gff_detect_width/;
use Moose;

with 'GBUtil::InputFile';

Readonly my @types =>  qw{methyl coverage};

sub BUILD{
    my ($self) = @_;
}

has wig_dir => (
    is => 'ro',
    required => 1,
);

has window => (
    is => 'ro', 
    lazy_build => 1,
    init_arg => 0,
);

sub _build_window{
    my $self = shift;
    return gff_detect_width($self->file);
}

has ctscore => (
    is => 'ro',
    default => 0,
    required => 0,
);

has [qw/_wig_files _meta_files/] => (
    is => 'rw',
    default => sub{ {} },
    init_arg => undef,
);

sub get_wig_file{
    my ($self, $feature, $type, $seq) = @_;
    my $f = $self->_wig_files()->{$feature}{$type}{$seq};
    if ($f){
        return $f;
    }
    return; # so that () in list context
}

sub set_wig_file{
    my ($self, $feature, $type, $seq, $filename) = @_;
    $self->_wig_files()->{$feature}{$type}{$seq} = $filename;
}

sub get_meta_file{
    my ($self, $feature, $type) = @_;
    my $f = $self->_meta_files()->{$feature}{$type};
    if ($f){
        return $f;
    }
    return; # so that () in list context
}

sub set_meta_file{
    my ($self, $feature, $type, $filename) = @_;
    $self->_meta_files()->{$feature}{$type} = $filename;
}

sub tracks{
    my ($self) = @_;
    my $h = $self->_meta_files;
    my @rv;
    for my $feature (sort keys %$h) {
        for my $type (@types) {
            push @rv, {
                meta    => $self->get_meta_file($feature, $type),
                feature => $feature,
                type    => $type,
            };
        }
    }
    return @rv;
}

sub upload_files{
    my ($self) = @_;
    return map { $_->{meta} } $self->tracks;
}

sub meta_files{
    my ($self) = @_;
    my $h = $self->_meta_files;
    my @rv;
    for my $feature (sort keys %$h) {
        for my $type (@types) {
            push @rv, $self->get_meta_file($feature, $type),
        }
    }
    return @rv;
}

has _features => (
    traits    => ['Hash'],
    is        => 'ro',
    isa       => 'HashRef[Str]',
    default   => sub { {} },
    handles   => {
        set_feature => 'set',
        has_feature => 'exists',
        features    => 'keys'
    },
);

has _sequences => (
    traits    => ['Hash'],
    is        => 'ro',
    isa       => 'HashRef[Str]',
    default   => sub { {} },
    handles   => {
        set_sequence => 'set',
        has_sequence => 'exists',
        sequences    => 'keys'
    },
);

#######################################################################
# methylation gff

# compile wiggle to base_directory and write write meta file
sub compile_wiggle{
    my $self = shift;
    my %opt = validate(@_, {
            feature        => 1, 
            track          => 1, # not sure what this is
            meta_file      => 1,
            files => {
                type => ARRAYREF,
            },
        });
    lock_keys(%opt);
    my $source = $self->source;

    my $loader = Bio::Graphics::Wiggle::Loader->new($self->wig_dir)
        or die "could not create loader";

    $loader->{trackname} = $opt{track} if defined $opt{track};

    for my $file (@{$opt{files}}) {
        say STDERR "compile_wiggle $opt{feature} $source " . basename($self->file);
        my $fh = IO::File->new($file) or die "could not open $file: $!";
        $loader->load($fh);
    }

    say STDERR "writing meta_file for $opt{feature} $source $opt{meta_file}";
    open my $metafh, '>', $opt{meta_file};
    print {$metafh} $loader->featurefile('gff3',$opt{feature},$source);
    close $metafh;
}

# returns
#   [ 
#     { feature => 'CG', type => 'methyl | coverage', meta => 'meta.gff },
#     ...
#   ]
sub convert{
    my $self = shift;

    make_path $self->staging_dir;

    my $window_size = $self->window;
    my %window_size_mismatch;

    my %wig_handles;

    my $p = GFF::Parser->new(file => $self->file, normalize => 0);

    my $ctscore = $self->ctscore;

    # convert gff to wig 
    while (defined(my $gff = $p->next())){
        my ($seq, $start, $feature, $end, $score, $c, $t) 
        = (lc($gff->sequence()), $gff->start(), $gff->feature(), $gff->end(), $gff->score(), $gff->get_column('c') // 0, $gff->get_column('t') // 0,);

        my $width = $end - $start + 1;

        # okay for a single window to be different size b/c the last one may be smaller
        if ($width != $window_size and $window_size_mismatch{$seq}++ > 1){
            croak "uneven widths. line $.\n$gff"; 
        }

        $score //= "0";
        if ($ctscore and defined($c) and defined($t)){ $score = ($c + $t) == 0 ? 0 : $c / ($c + $t); }

        # create wigs if ! exist and write header
        for my $type (@types) {
            if (!$self->get_meta_file($feature, $type)){
                $self->set_meta_file($feature, $type, 
                    $self->make_meta_filename(
                        feature    => $feature,
                        type       => $type
                    )
                );
            }
            if (!$self->get_wig_file($feature, $type, $seq)){
                my $filename = $self->make_wig_filename(
                    seq        => $seq,
                    feature    => $feature,
                    type       => $type
                );
                open my $fh, '>', $filename;

                $self->set_wig_file($feature, $type, $seq, $filename);
                if (! $self->has_feature($feature)){
                    $self->set_feature($feature,1);
                }
                if (! $self->has_sequence($seq)){
                    $self->set_sequence($seq,1);
                }

                $wig_handles{$feature}{$type}{$seq} = $fh;

                say STDERR "created " .  $self->source()  . " $seq $feature $type w$window_size";

                say $fh qq{variableStep chrom=$seq  span=$window_size};
                # say $fh qq{variableStep name="$feature-$type" chrom=$seq  span=$detected_width};
            }
        }

        say {$wig_handles{$feature}{methyl}{$seq}}   sprintf("%d\t%4f", $start, $score);
        say {$wig_handles{$feature}{coverage}{$seq}} "$start\t", $c + $t;
    }

    # close all handles
    for my $f (keys %wig_handles) {
        for my $type (qw/methyl coverage/) {
            for my $s (keys %{$wig_handles{$f}{$type}}) {
                close $wig_handles{$f}{$type}{$s}
            }
        }
    }

    # compile with compile_wiggle 
    for my $feature ($self->features){
        for my $type (qw/methyl coverage/) {
            say STDERR "done reading " . $self->title;

            $self->compile_wiggle(
                files     => [map { $self->get_wig_file($feature, $type, $_) } $self->sequences],
                track     => $self->source() . "-$feature-$type-$window_size",
                feature   => "$feature-$type-$window_size",
                meta_file => $self->get_meta_file($feature, $type),
            );
        }
    }
}

#######################################################################

sub _make_tmpname{
    my $self = shift;
    my $abs = rel2abs($self->file);
    return join "_", grep { $_ ne '' } split /\//, $abs;
}

#$meta_files{$feature}{$type} = sprintf($meta_filename_template, $feature, $type);
sub make_meta_filename{
    my $self = shift;
    my %opt = validate(@_, { feature => 1, type => 1, });
    lock_keys(%opt);
    return sprintf(
        catfile($self->staging_dir, $self->_make_tmpname($self->file)) . "-%s-%s.meta.gff",
        $opt{feature}, $opt{type}
    );
}

# my $filename = sprintf($wig_filename_template, $seq, $feature, $type);
sub make_wig_filename{
    my $self = shift;
    my %opt = validate(@_, { feature => 1, type => 1, seq => 1, });
    return sprintf(
        catfile($self->staging_dir, $self->_make_tmpname($self->file)) . "-%s-%s.%s.wig",
        $opt{seq}, $opt{feature}, $opt{type}
    );
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
