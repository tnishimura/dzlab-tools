#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use feature 'say';
use Archive::Tar;
use Cwd qw/getcwd/;
use Data::Dumper;
use File::Basename qw/basename dirname/;
use File::Copy;
use File::Find;
use File::Path qw/make_path remove_tree/;
use File::Spec::Functions qw/rel2abs canonpath catdir catfile updir/;
use File::Temp qw/tempdir/;
use Getopt::Long;
use Git::Repository;
use Pod::Usage;
use POSIX qw/strftime/;
use Scalar::Util qw/looks_like_number/;
use Template;
use YAML qw/Load DumpFile LoadFile/;

$0 = basename($0);

my $result = GetOptions (
    "config|c=s" => \(my $config_file = catfile($ENV{HOME}, "etc", "dzlab-tools.conf")),
    "debug|d" => \my $debug,
);
if (!$result || ! -f $config_file){
    say "usage: $0 [--config dzlab-tools.conf]";
    exit 1;
}

#######################################################################
# read $HOME/etc/dzlab-tools.conf

my $config = LoadFile($config_file);
my          ($makensis, $masterdir, $earliest, $releasedir) = 
@{$config}{qw/makensis   masterdir   earliest   releasedir/};

$makensis //= "/usr/bin/makensis";
die "masterdir in $config_file not a directory?" if (! -d $masterdir);
die "releasedir in $config_file not a directory?" if (! -d $releasedir);
die "earliest in $config_file not a number?" if (! looks_like_number($earliest));

#######################################################################
# read DATA section

my ($ignore_dirs, $installer_nsi_template, $dzlab_check_pl_template) = 
@{Load(do {local $/; scalar <DATA>})}{qw/ignore_dirs installer_nsi_template dzlab_check_pl_template/};

#######################################################################

{
    my $masterrepo = Git::Repository->new( git_dir => $masterdir );

    sub _all_tags{
        grep { /^\d+\.\d+\.\d+$/ } $masterrepo->run('tag');
    }

    sub all_recent_tags{
        grep { 
            my $d = tagdate($_);
            looks_like_number($d) && $d > $earliest
        } _all_tags;
    }

    sub tagdate{
        my $tag = shift;
        my $timestamp = $masterrepo->run(show => qw/-s --format=%at/, $tag);
        chomp $timestamp;
        return $timestamp;
    }
}

sub clone_tag_to_dir{
    my ($tag, $dir) = @_;

    Git::Repository->run(clone => $masterdir, $dir);

    my $localrepo = Git::Repository->new(work_tree => $dir);
    $localrepo->run(checkout => $tag);
}

sub should_ignore_dir{
    my $dir = shift;
    return 0 < grep { $dir eq $_ } @$ignore_dirs;
}

sub should_ignore_file{
    my $file = shift;
    return 1 if $file =~ /^\./ || $file eq $0 || $file =~ /\.in$/;
}

sub dump_version_file{
    my ($file, $tag) = @_;
    DumpFile($file, {
            buildtime => strftime("%Y/%m/%d %H:%M:%S %Z", localtime(time())),
            committime => strftime("%Y/%m/%d %H:%M:%S %Z", localtime(tagdate($tag))),
            version => $tag,
        });
}

sub build_linux{
    my ($tag, $dest) = @_;

    my $tempdir  = tempdir(CLEANUP => !$debug);
    my $builddir = catdir($tempdir, 'dzlab-tools');
    make_path($builddir);

    clone_tag_to_dir($tag,$builddir);
    dump_version_file catfile($builddir, "VERSION"), $tag;

    my $cwd = getcwd();

    chdir $tempdir; say "*changed to $tempdir";

    my @files;
    my $tar = Archive::Tar->new;

    find( sub {
            # $File::Find::name - filename relative to pwd
            # $File::Find::dir  - dirname relative to pwd, eq getcwd()
            # $_                - filename relative to $File::Find::dir

            if (-d && should_ignore_dir $_){
                $File::Find::prune = 1;
            }
            elsif (-f && ! should_ignore_file $_){
                # push to array first, b/c find() chdirs to each dir by default
                push @files, canonpath $File::Find::name;
            }
        }, '.');

    $tar->add_files(@files);
    $tar->write($dest, COMPRESS_GZIP);

    chdir $cwd;
}

sub build_windows{
    my ($tag, $dest, $basename) = @_;

    my $tempdir  = tempdir(CLEANUP => !$debug);
    my $tt = Template->new() || die "$Template::ERROR\n";

    clone_tag_to_dir($tag,$tempdir);
    dump_version_file catfile($tempdir, "VERSION"), $tag;

    my $cwd = getcwd();
    chdir $tempdir; 

    $tt->process(\$installer_nsi_template, {version => $tag}, "installer.nsi")
    || die $tt->error(), "\n";

    $tt->process(\$dzlab_check_pl_template, {version => $tag}, "dzlab-check.pl")
    || die $tt->error(), "\n";

    if (0 == system("$makensis installer.nsi")){
        move $basename, $dest;
    }
    else{
        die "couldn't execute makensis?";
    }

    chdir $cwd;
}

for my $tag (all_recent_tags) {
    my $windows_name = "dzlab-tools-$tag.exe";
    my $linux_name   = "dzlab-tools-$tag.linux.tar.gz";

    my $windows_dest = catfile($releasedir, $windows_name);
    my $linux_dest   = catfile($releasedir, $linux_name);

    # windows 
    if (! -e $windows_dest && -x $makensis){
        build_windows($tag, $windows_dest, $windows_name);
    }

    # linux
    if (! -f $linux_dest){
        build_linux($tag, $linux_dest);
    }
}

__DATA__
---
ignore_dirs:
  - util
  - utilities
  - .git
  - t
  - tmp

dzlab_check_pl_template: |-
  #!/usr/bin/env perl
  use strict;
  use warnings;
  
  print "DZLab-Tools [% version %] installed! (April 11, 2011)\n";
  exit 0;

installer_nsi_template: |-
  outFile "dzlab-tools-[% version %].exe"
   
  InstallDir C:\dzlab-tools
  !define UNINSTALLER $INSTDIR\uninstaller.exe
  
  # execute unintaller if it already exists
  # http://nsis.sourceforge.net/When_I_use_ExecWait_uninstaller.exe_it_doesn't_wait_for_the_uninstaller
  Function .onInit
  ExecWait '"${UNINSTALLER}" _?=$INSTDIR'
  FunctionEnd
   
  section
  setOutPath $INSTDIR
  file /r /x utilities /x .* /x .git /x *.exe /x *.nsi /x t /x tmp /x *.gff /x _build* *
  writeUninstaller ${UNINSTALLER}
  CreateDirectory "$SMPROGRAMS\Zilberman Lab"
  CreateShortCut "$SMPROGRAMS\Zilberman Lab\Help.lnk" "$INSTDIR\help\index.html"
  sectionEnd
   
  section "Uninstall"
  SetAutoClose true
  Delete "$SMPROGRAMS\Zilberman Lab\Help.lnk"
  delete ${UNINSTALLER}
  rmdir /r $INSTDIR\*
  sectionEnd
