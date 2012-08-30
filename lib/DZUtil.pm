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
use File::Temp qw/tempdir tempfile/;
use File::Path qw/make_path/;
use List::MoreUtils qw/all/;
use POSIX qw/strftime/;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::Uncompress::Bunzip2 qw(bunzip2 $Bunzip2Error) ;
use Scalar::Util qw/looks_like_number/;
use List::Util qw//;
use List::MoreUtils qw//;
use YAML qw/LoadFile DumpFile/;

require Exporter;
our @ISA = qw(Exporter);

our @EXPORT_OK = qw(localize reverse_complement common_suffix common_prefix
mfor basename_prefix fastq_read_length timestamp datestamp overlap chext
split_names open_maybe_compressed fastq_convert_read_header c2t g2a rc_c2t rc_g2a
numdiff safediv safemethyl clean_basename open_cached close_cached_all downsample
approximate_line_count memofile memodo gimmetmpdir split_file combine_csv);
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

sub strip_common_prefix{
    my $pre = common_prefix(@_);
    return map { s/^\Q$pre\E// } @_;
}

sub strip_common_suffix{
    my $suf = common_suffix(@_);
    return map { s/\Q$suf\E$// } @_;
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

# Code  Bases  RC  C2T  G2A
# A     A      T   A    A
# B     CGT    V   K    H
# C     C      G   T    C
# D     AGT    H   D    W
# G     G      C   G    A
# H     ACT    D   W    H
# K     GT     M   K    W
# M     AC     K   W    M
# N     ACGT   N   D    H
# R     AG     Y   R    A
# S     CG     S   K    M
# T     T      A   T    T
# V     ACG    B   D    M
# W     AT     W   W    W
# Y     CT     R   T    Y
# rc: tr/abcdghkmrtvyABCDGHKMRTVY/tvghcdmkyabrTVGHCDMKYABR/
# c2t: tr/bchmnsvyBCHMNSVY/ktwwdkdtKTWWDKDT/
# g2a: tr/bdgknrsvBDGKNRSV/hwawhammHWAWHAMM/

sub reverse_complement{
    my $string = shift;

    $string = reverse $string;
    $string =~ tr/abcdghkmrtvyABCDGHKMRTVY/tvghcdmkyabrTVGHCDMKYABR/;
    return $string;
}

sub c2t{
    my $string = shift;
    $string =~ tr/bchmnsvyBCHMNSVY/ktwwdkdtKTWWDKDT/;
    return $string;
}

sub rc_c2t{
    my $string = shift;
    return c2t reverse_complement $string;
}

sub g2a{
    my $string = shift;
    $string =~ tr/bdgknrsvBDGKNRSV/hwawhammHWAWHAMM/;
    return $string;
}

sub rc_g2a{
    my $string = shift;
    return g2a reverse_complement $string;
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

sub downsample{
    my ($file, $percent, $multiple, $tempdir) = @_;
    $multiple //= 1;
    if (defined $percent && looks_like_number $percent && $percent < 0 || $percent > 1){
        croak "downsample: \$percent needs to be 0<=p<=1";
    }

    $tempdir //= tempdir(CLEANUP => 1);
    my $tempfile = catfile($tempdir, basename($file) . ".downsample$percent");

    my $infh = open_maybe_compressed $file;
    open my $outfh, '>', $tempfile;

    my $sampled = 0;
    my $total = 0;
    while (defined(my $line = <$infh>)){
        my @others = map { (eof $infh && last) || scalar <$infh> } (2 .. $multiple);
        if (rand() < $percent){
            for my $l ($line, @others) {
                print $outfh $l;
            }
            ++$sampled;
        }
        ++$total;
    }
    close $outfh;
    close $infh;
    if (wantarray){
        return $tempfile, $sampled * $multiple, $total * $multiple;
    }
    else{
        return $tempfile;
    }
}


sub approximate_line_count{
    my ($file, $max) = @_;
    my $size = [stat($file)]->[7];
    $max //= 100_000;

    my @lengths;
    my $count = 0;
    open my $in, '<', $file;
    while (defined(my $line = <$in>) && ++$count <= $max){
        chomp $line;
        push @lengths, length $line;
    }
    close $in;
    my $mean = List::Util::sum(@lengths)/scalar(@lengths);

    return int $size/$mean;
}

# flatten a path and return it.  if there's a second argument, catfile() onto it.
# memofile("/home/kitten/meow.txt", "/tmp") => /tmp/home,kitten,meow.txt
sub memofile{
    my ($path, $tmpdir) = @_;
    my (undef,$directories,$file) = File::Spec::Functions::splitpath(File::Spec::Functions::rel2abs($path)); 

    my @directory_parts = File::Spec::Functions::splitdir($directories);
    @directory_parts = grep { $_ ne '' } @directory_parts; 

    return catfile($tmpdir, join ",", @directory_parts, $file);
}

# do $code if memofile $file doesn't exist (or forced)
sub memodo{
    my ($file, $code, $force) = @_;
    if (! -f $file || $force){
        my $results = $code->();
        DumpFile($file, $results);
        return $results;
    }
    else{
        return LoadFile($file);
    }
}

# without an argument, create and return a tmpdir.
# with an argument, make sure it exists and that its a directory.
sub gimmetmpdir{
    my ($tmpdir) = @_;
    $tmpdir //= tempdir(CLEANUP => 1);

    if (-e $tmpdir){
        if (! -d $tmpdir){
            die "$tmpdir not a dir";
        }
    }
    else{
        make_path $tmpdir;
        if (! -d $tmpdir){
            die "can't make $tmpdir";
        }
    }
    return $tmpdir;
}

use Params::Validate qw/:all/;
sub split_file{
    my %opt = validate(@_, {
            file => {
                type => SCALAR,
                callbacks => {
                    'file must exist' => sub { -f shift },
                },
                optional => 1,
            }, 
            outdir   => { 
                type => SCALAR, callbacks => {
                    'directory must exist' => sub { -d shift },
                },
                optional => 1,
            },
            size     => { type => SCALAR, default => 2**30, }, # gigabyte 
            multiple => { type => SCALAR, default => 1, },
            skip     => { type => SCALAR, default => 0, }
        });

    my $prefix = $opt{outdir} ? catfile($opt{outdir}, basename($opt{file})) : basename($opt{file});
    $prefix .= ".prefix";

    my $hashfile = $prefix . ".MD5SUM";

    if (-f $hashfile){
        my $md5sums = LoadFile($hashfile);
        my $ok = 1;
        while (my ($file,$md5) = each %$md5sums) {
            if (md5sum($file) ne $md5){
                $ok = 0;
            }
        }
        if ($ok == 1){
            return keys %$md5sums;
        }
    }


    my $suffix = "aaa";
    my $bytes_written = 0;;

    my @split;

    open my $reader, '<', $opt{file};
    open my $writer, '>', "$prefix.$suffix";
    push @split, "$prefix.$suffix";

    while (! eof $reader){
        my @lines = map { scalar readline $reader } (1 .. $opt{multiple});
        my $to_write = List::Util::sum(map { defined() ? length($_) : 0 } @lines);

        if (List::MoreUtils::notall {defined $_} @lines){
            carp "$opt{file} was not a multiple of $opt{multiple}";
        }

        if ($opt{skip}){
            for (1 .. $opt{skip}) {
                scalar readline $reader;
            }
        }

        for my $line (@lines) {
            print $writer $line;
        }

        $bytes_written += $to_write;

        if ($bytes_written > $opt{size}){
            close $writer;
            $suffix++;
            open $writer, '>', "$prefix.$suffix";
            push @split, "$prefix.$suffix";
            $bytes_written = 0;
        }
    }

    close $writer;
    close $reader;

    @split = map { File::Spec::Functions::rel2abs($_) } @split;
    DumpFile($hashfile, { map { $_ => md5sum($_) } @split });

    return @split;

}

# combine multiple CSV files so that only the first one has a header.
# croak when the # of columns mismatch.
sub combine_csv{
    my ($outfile, @infiles) = @_;
    my $first_file = 1;

    my $colcount;

    my $tmpout = $outfile . ".tmp";

    open my $outfh, '>:crlf', $tmpout;

    for my $inf (@infiles) {
        open my $infh, '<:crlf', $inf;
        my $first_line = 1;

        while (defined(my $line = <$infh>)){
            chomp $line;

            # check that it has correct number of columns
            if ($first_line){
                my $file_colcount = scalar(split /\t/, $line);
                $colcount = $colcount // $file_colcount;
                if ($colcount != $file_colcount){
                    croak "$inf colcount mismatch?";
                }

                # and print it if its also the first file
                if ($first_file){
                    say {$outfh} $line;
                }
            }
            else{
                # if its not first line, always print
                say {$outfh} $line;
            }

            $first_line = 0;
            $first_file = 0;
        }

        close $infh;
    }

    close $outfh;

    rename $tmpout, $outfile;
}

1;

