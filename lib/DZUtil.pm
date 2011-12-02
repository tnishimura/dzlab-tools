package DZUtil;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;
use File::Spec::Functions;
use File::Basename;
use List::MoreUtils qw/all/;
#use Config::General qw(ParseConfig);
use POSIX qw/strftime/;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::Uncompress::Bunzip2 qw(bunzip2 $Bunzip2Error) ;

require Exporter;
our @ISA = qw(Exporter);

our @EXPORT_OK = qw(localize reverse_complement common_suffix common_prefix
mfor basename_prefix fastq_read_length timestamp datestamp overlap chext
split_names base_match open_maybe_compressed fastq_convert_read_header c2t
numdiff safediv safemethyl clean_basename open_cached close_cached_all);
our @EXPORT = qw();

sub clean_basename{
    my ($path, @exts) = @_;
    defined $path or croak "need filename";
    $path = basename($path, @exts);
    $path =~ s/[^\s.\w-]//g; # get rid of funny chars
    $path =~ s/\s+/_/g;     # spaces to underscores
    return $path;
}

{
    my %cached_filehandles;
    sub open_cached{
        my ($mode, $file) = @_;
        if (exists $cached_filehandles{$file}){
            return $cached_filehandles{$file};
        }
        else{
            open my($fh), $mode, $file;
            $cached_filehandles{$file} = $fh;
            return $fh;
        }
    }
    sub close_cached_all{
        close $_ for values %cached_filehandles;
        %cached_filehandles = ();
    }
}

sub c2t{
    (my $c2t = $_[0]) =~ s/C/T/gi;
    return $c2t;
}

sub numdiff{
    my ($x,$y)=@_;
    my @x_split = split //, $x;
    my @y_split = split //, $y;
    croak "len mismatch" unless @x_split == @y_split;
    my $total = 0;
    for (0..$#x_split){
        if ($x_split[$_] ne $y_split[$_]){
            $total += 1;
        }
    }
    return $total;
}


=head2 chext("/etc/test.txt", "newext")

return filename with new extension

=cut

sub chext{
    my ($fullpath, $ext) = @_;
    my ($filename, $path) = fileparse($fullpath);
    $filename =~ s/\.[^.]+$//;
    return catfile($path,$filename) . ".$ext";
}

=head2 split_names("/home/test/path.txt", qw/chr1 chr2 chr3/)

Yields (/home/test/path-chr1.txt /home/test/path-chr2.txt /home/test/path-chr2.txt)

=cut

sub split_names{
    my ($fullpath, @splits) = @_;
    my ($filename, $path) = fileparse($fullpath);
    my $ext = ($filename =~ s/\.([^.]+)$//) ? $1 : "";

    return map { catfile($path, $filename) . "-$_" . ($ext ? ".$ext" : "") } @splits;
}

=head2 overlap([start1,end1], [start2,end2])

Return overlap between two ranges, or 0 if not overlapping.

=cut

use List::Util qw/max min first/;

sub overlap{
    my ($x,$y) = @_;
    my $start1 = min($x->[0], $x->[1]);
    my $end1   = max($x->[0], $x->[1]);

    my $start2 = min($y->[0], $y->[1]);
    my $end2   = max($y->[0], $y->[1]);

    if ($end1 >= $start2 && $end2 >= $start1){
        return min($end1, $end2) - max($start1, $start2)  + 1;
    }
    else {
        return 0;
    }
}


my %iupac = (
    G => [qw/ G /],
    A => [qw/ A /],
    T => [qw/ T /],
    C => [qw/ C /],
    R => [qw/ G A /],
    Y => [qw/ T C /],
    M => [qw/ A C /],
    K => [qw/ G T /],
    S => [qw/ G C /],
    W => [qw/ A T /],
    H => [qw/ A C T Y M W /],
    B => [qw/ G T C Y K S /],
    V => [qw/ G C A R M S /],
    D => [qw/ G A T R K W /],
    N => [qw/ G A T C R Y M K S W H B V D N/],
);

# base_match - return true if base1 is compatible with base2. 

sub base_match{
    my ($base1, $base2, $both) = @_;

    $base1 = uc $base1;
    $base2 = uc $base2;

    return first{ $base1 eq $_ } (@{$iupac{$base2}});
}

sub timestamp{ return strftime("%Y%m%d-%H%M",localtime); }
sub datestamp{ return strftime("%Y%m%d",localtime); }

sub fastq_read_length{
    my $prefix = shift;
    if (-f $prefix){
        return _fastq_read_length_single($prefix);
    }
    # assume a prefix
    else{
        my @sizes = map { 
            _fastq_read_length_single($_) 
        } glob(join " ", catfile($prefix, "*.gz"), catfile($prefix, "*.bz2"));
        if (@sizes == 0){
            croak "no files in $prefix";
        }
        else{
            my $first = shift @sizes;
            if (all { $first == $_ } @sizes){
                return $first;
            }
            else{
                croak "uneven read lengths in $prefix";
            }
        }
    }
}

sub _fastq_read_length_single{
    my $filename = shift;
    my $fh = open_maybe_compressed($filename);
    <$fh>;
    my $line = <$fh>;
    close $fh;

    if (defined $line){
        #chomp $line;
        $line =~ tr/\n\r//d;
        return length $line;
    }
    return;
}

sub fastq_convert_read_header{
    my $line = $_[0];
    $line =~ s/^@//;
    my ($first, $second) = split ' ', $line;

    # most common
    if ($first =~ m{#0/[12]$}){
        return $first;
    }
    # new solexa format, 6/20/2011
    elsif (defined $second && $second =~ /^([12])/){
        return "$first#/$1";
    }
    # handle older reads with no read id's 
    else{
        return "$first#/1";
    }
}

=head2 open_maybe_compressed

return a file handle for a file which may be compressed

=cut 
sub open_maybe_compressed{
    my $filename = shift;
    if ($filename =~ /\.gz$/i){
        return new IO::Uncompress::Gunzip $filename 
            or croak "IO::Uncompress::Gunzip failed: $GunzipError\n";
    }
    elsif ($filename =~ /\.bz2$/i){
        return new IO::Uncompress::Bunzip2 $filename 
            or croak "IO::Uncompress::Bunzip2 failed: $Bunzip2Error\n";
    }
    else {
        open my $fh, '<', $filename;
        return $fh;
    }
}

#sub read_conf{
#    my ($index) = grep { '--conf' eq $main::ARGV[$_] } (0 .. $#main::ARGV);
#    return () if ! defined $index;
#
#    my (undef, $file) = splice @main::ARGV, $index, 2;
#    return () if ! $file;
#
#    croak "--conf $file not a file?" if ! -f $file;
#
#    return ParseConfig("$file");
#}

sub basename_prefix{
    my $file = shift;
    my $base = basename($file);
    if ($base =~ /([^.]*?)\./){
        return $1;
    }
    return $base;
}


sub common_prefix{
    my $min = min map { length $_ } @_;
    my @arrayed = map { [split //, $_] } @_;
    
    my @accum;
    for my $i (0 .. $min-1) {
        my @ith_chars = map { $_->[$i] } @arrayed;
        my $first = pop @ith_chars;
        if (grep { $first ne $_ } @ith_chars){
            last;
        } 
        else {
            push @accum, $first;
        }
    }
    return join "", @accum;
}

sub common_suffix{
    return scalar reverse common_prefix map { scalar reverse $_ } @_;
}

sub mfor {
    my $code = pop;
    my @arrays = @_;
    my $minlength = min map { scalar @$_ } @arrays;

    for my $i (0 .. $minlength - 1) {
        $code->(map { $_->[$i] } @arrays);
    }
}


sub reverse_complement{
    my $string = shift;

    $string = reverse $string;
    $string =~ tr/acgtACGT/tgcaTGCA/;
    return $string;
}

use LWP::Simple;
use File::Temp qw/tempdir mktemp/;

sub localize{
    my ($file_or_url, $dir, $overwrite) = @_;
    if (!defined $dir || ! -d $dir){
        $dir = tempdir(CLEANUP => 1);
    }

    if ($file_or_url =~ m{^(http|ftp)://}){
        my $storefile = catfile($dir,basename($file_or_url));
        my $tmpfile = mktemp($storefile . "XXXX");
        if (-f $storefile){
            if ($overwrite){
                carp "localize(): $storefile already exists? overwriting";
            }
            else{
                #carp "localize(): $storefile already exists? NOT overwriting";
                return $storefile;
            }
        }

        if (200 == mirror($file_or_url, $tmpfile)){
            rename $tmpfile, $storefile;
            return $storefile;
        }
        else {
            croak "couldn't grab $file_or_url into $storefile";
        }
    }
    else {
        return $file_or_url;
    }
}

sub safediv{
    my ($num, $den) = @_;
    if (defined $den && $den > 0){
        return $num/$den;
    }
    return 0;
}

sub safemethyl{
    my ($c, $t) = @_;
    if ($c + $t > 0){
        return $c / ($c + $t);
    }
    else{
        return 0;
    }
}

1;


