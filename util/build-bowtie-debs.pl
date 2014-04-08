#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Cwd qw/getcwd/;
use Data::Dumper;
use Digest::MD5;
use File::Copy;
use File::Path qw/make_path remove_tree/;
use File::Temp qw/mkdtemp/;
use LWP::Simple;
use YAML qw/Load/;
use Getopt::Long;
use Pod::Usage;

my $result = GetOptions (
    "debug|d" => \(my $debug),
);
pod2usage(-verbose => 2, -noperldoc => 1) if ! $result;

if (0 != system("which fpm")){
    say "can't find fpm in \$PATH?";
    exit 1;
}

my $packages = Load(do {local $/; scalar <DATA>});

make_path("build");
chdir("build");

for my $package (@$packages) {
    my $name    = $package->{name};
    my $version = $package->{version};
    my $url     = $package->{url};
    my $license = $package->{license};
    my $tarball = $package->{tarball};
    my $md5sum  = $package->{md5sum};
    my $bin_files   = $package->{bin_files};

    if (! -f $tarball || md5sum($tarball) ne $md5sum){
        getstore($url, $tarball);
    }
    if (! -f $tarball || md5sum($tarball) ne $md5sum){
        die "couldn't get $tarball?";
    }

    if ($tarball =~ /\.zip$/){
        system("unzip -q -o $tarball");
    }
    elsif ($tarball =~ /\.tar.gz$/){
        system("tar xzf $tarball");
    }

    my $tmpdir = mkdtemp("build-$name-XXXXXX");
    my $bin = "$tmpdir/usr/bin/";
    my $share = "$tmpdir/usr/share/doc/$name/";
    make_path($bin);
    make_path($share);

    for my $f (@$bin_files) {
        copy($f, $bin);
        chmod 0755, glob("$bin/*");
    }

    copy($license, $share);

    system(qq{fpm -a x86_64 -f --url http://dzlab.pmb.berkeley.edu --maintainer dzlab --vendor dzlab --description "For DZLab internal use only" -s dir -t deb -v $version -n $name -C $tmpdir .});
    remove_tree($tmpdir) unless $debug;
}

#######################################################################

sub md5sum{
    my $file = shift;
    # md5sum doesn't support addfile(filename)?
    open my $fh, '<', $file;
    binmode($fh);
    my $rv = Digest::MD5->new->addfile($fh)->hexdigest();
    close $fh;
    return $rv;
}

__DATA__
---

- name : bowtie
  version : 0.12.9
  url : http://downloads.sourceforge.net/project/bowtie-bio/bowtie/0.12.9/bowtie-0.12.9-linux-x86_64.zip
  tarball : bowtie-0.12.9-linux-x86_64.zip
  md5sum : 3489b9ca2e398011f47baa75f9f45122
  license: bowtie-0.12.9/COPYING
  bin_files : 
      - ./bowtie-0.12.9/bowtie
      - ./bowtie-0.12.9/bowtie-build
      - ./bowtie-0.12.9/bowtie-build-debug
      - ./bowtie-0.12.9/bowtie-debug
      - ./bowtie-0.12.9/bowtie-inspect
      - ./bowtie-0.12.9/bowtie-inspect-debug

- name : bowtie
  version : 1.0.0
  url : http://downloads.sourceforge.net/project/bowtie-bio/bowtie/1.0.0/bowtie-1.0.0-linux-x86_64.zip
  tarball : bowtie-1.0.0-linux-x86_64.zip
  md5sum : c1c3b6599fd0d92213b130cbccfc5e92
  license: bowtie-1.0.0/LICENSE
  bin_files : 
      - ./bowtie-1.0.0/bowtie
      - ./bowtie-1.0.0/bowtie-build
      - ./bowtie-1.0.0/bowtie-build-debug
      - ./bowtie-1.0.0/bowtie-debug
      - ./bowtie-1.0.0/bowtie-inspect
      - ./bowtie-1.0.0/bowtie-inspect-debug

- name : bowtie2
  version : 2.1.0
  url : http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zip
  tarball : bowtie2-2.1.0-linux-x86_64.zip
  md5sum : c00c168f8825382896cee967e191d71b
  license: bowtie2-2.1.0/LICENSE
  bin_files : 
      - ./bowtie2-2.1.0/bowtie2
      - ./bowtie2-2.1.0/bowtie2-align
      - ./bowtie2-2.1.0/bowtie2-align-debug
      - ./bowtie2-2.1.0/bowtie2-build
      - ./bowtie2-2.1.0/bowtie2-build-debug
      - ./bowtie2-2.1.0/bowtie2-inspect
      - ./bowtie2-2.1.0/bowtie2-inspect-debug

- name : cufflinks
  version : 1.3.0
  url : http://cufflinks.cbcb.umd.edu/downloads/cufflinks-1.3.0.Linux_x86_64.tar.gz
  tarball : cufflinks-1.3.0.Linux_x86_64.tar.gz
  md5sum : 12e6b955aadca29e02570ddea83ad9b6
  license: cufflinks-1.3.0.Linux_x86_64/LICENSE
  bin_files : 
      - ./cufflinks-1.3.0.Linux_x86_64/cuffcompare
      - ./cufflinks-1.3.0.Linux_x86_64/cuffdiff
      - ./cufflinks-1.3.0.Linux_x86_64/cufflinks
      - ./cufflinks-1.3.0.Linux_x86_64/cuffmerge
      - ./cufflinks-1.3.0.Linux_x86_64/gffread
      - ./cufflinks-1.3.0.Linux_x86_64/gtf_to_sam

- name : cufflinks
  version : 2.2.0
  url : http://cufflinks.cbcb.umd.edu/downloads/cufflinks-2.2.0.Linux_x86_64.tar.gz
  tarball : cufflinks-2.2.0.Linux_x86_64.tar.gz
  md5sum : 589b8b9b4dbb036719bb75682fced057
  license: cufflinks-2.2.0.Linux_x86_64/LICENSE
  bin_files : 
      - ./cufflinks-2.2.0.Linux_x86_64/cuffcompare
      - ./cufflinks-2.2.0.Linux_x86_64/cuffdiff
      - ./cufflinks-2.2.0.Linux_x86_64/cufflinks
      - ./cufflinks-2.2.0.Linux_x86_64/cuffmerge
      - ./cufflinks-2.2.0.Linux_x86_64/gffread
      - ./cufflinks-2.2.0.Linux_x86_64/gtf_to_sam

- name : tophat
  version : 1.4.1
  url : http://tophat.cbcb.umd.edu/downloads/tophat-1.4.1.Linux_x86_64.tar.gz
  tarball : tophat-1.4.1.Linux_x86_64.tar.gz
  md5sum    : f4a3243551c8a5a19b7a1ff115c18447
  license: tophat-1.4.1.Linux_x86_64/COPYING
  bin_files : 
      - ./tophat-1.4.1.Linux_x86_64/bam2fastx
      - ./tophat-1.4.1.Linux_x86_64/bam_merge
      - ./tophat-1.4.1.Linux_x86_64/bed_to_juncs
      - ./tophat-1.4.1.Linux_x86_64/closure_juncs
      - ./tophat-1.4.1.Linux_x86_64/contig_to_chr_coords
      - ./tophat-1.4.1.Linux_x86_64/extract_reads
      - ./tophat-1.4.1.Linux_x86_64/fix_map_ordering
      - ./tophat-1.4.1.Linux_x86_64/gtf_juncs
      - ./tophat-1.4.1.Linux_x86_64/gtf_to_fasta
      - ./tophat-1.4.1.Linux_x86_64/juncs_db
      - ./tophat-1.4.1.Linux_x86_64/library_stats
      - ./tophat-1.4.1.Linux_x86_64/long_spanning_reads
      - ./tophat-1.4.1.Linux_x86_64/map2gtf
      - ./tophat-1.4.1.Linux_x86_64/mask_sam
      - ./tophat-1.4.1.Linux_x86_64/prep_reads
      - ./tophat-1.4.1.Linux_x86_64/sam_juncs
      - ./tophat-1.4.1.Linux_x86_64/segment_juncs
      - ./tophat-1.4.1.Linux_x86_64/sra_to_solid
      - ./tophat-1.4.1.Linux_x86_64/tophat
      - ./tophat-1.4.1.Linux_x86_64/tophat_reports
      - ./tophat-1.4.1.Linux_x86_64/wiggles

- name : tophat
  version : 2.0.10
  url : http://tophat.cbcb.umd.edu/downloads/tophat-2.0.10.Linux_x86_64.tar.gz
  tarball : tophat-2.0.10.Linux_x86_64.tar.gz
  md5sum : af0a81b35af9a0c490e4ff748c22e182
  license: tophat-2.0.10.Linux_x86_64/COPYING
  bin_files : 
      - ./tophat-2.0.10.Linux_x86_64/bam2fastx
      - ./tophat-2.0.10.Linux_x86_64/bam_merge
      - ./tophat-2.0.10.Linux_x86_64/bed_to_juncs
      - ./tophat-2.0.10.Linux_x86_64/contig_to_chr_coords
      - ./tophat-2.0.10.Linux_x86_64/fix_map_ordering
      - ./tophat-2.0.10.Linux_x86_64/gtf_juncs
      - ./tophat-2.0.10.Linux_x86_64/gtf_to_fasta
      - ./tophat-2.0.10.Linux_x86_64/juncs_db
      - ./tophat-2.0.10.Linux_x86_64/long_spanning_reads
      - ./tophat-2.0.10.Linux_x86_64/map2gtf
      - ./tophat-2.0.10.Linux_x86_64/prep_reads
      - ./tophat-2.0.10.Linux_x86_64/sam_juncs
      - ./tophat-2.0.10.Linux_x86_64/segment_juncs
      - ./tophat-2.0.10.Linux_x86_64/sra_to_solid
      - ./tophat-2.0.10.Linux_x86_64/tophat
      - ./tophat-2.0.10.Linux_x86_64/tophat2
      - ./tophat-2.0.10.Linux_x86_64/tophat-fusion-post
      - ./tophat-2.0.10.Linux_x86_64/tophat_reports
