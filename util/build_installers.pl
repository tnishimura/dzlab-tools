#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use IO::All;
use Scalar::Util qw/looks_like_number/;
use FindBin;
use Pod::Usage;
use Getopt::Long;
use Cwd qw/getcwd/;
use File::Path qw/make_path remove_tree/;
use File::Spec::Functions qw/catdir catfile/;
use File::Copy;

my $result = GetOptions (
    "help"    => \my $help,
    "build-dir|b=s"   => \(my $build_dir   = "$ENV{HOME}/_build_dzlab_tools"),
    "release-dir|r=s" => \(my $release_dir = "$ENV{HOME}/_release_dzlab_tools"),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result);  

make_path($build_dir);
make_path($release_dir);

my $dzlab_tools_master     = "$FindBin::Bin/..";
my $installer_nsi_template = "$FindBin::Bin/installer.nsi.tt2";
my $dzlab_check_template   = "$FindBin::Bin/dzlab-check.pl..tt2";

die unless -f $installer_nsi_template;
chdir $dzlab_tools_master;

for my $tag (io("git tag |")->chomp->getlines()){
    if (tag_ok($tag)){
        # dzlab-tools-1.5.03.exe		 
        # dzlab-tools-1.5.03.linux.tar.gz  
        my $tag_build_dir              = catdir($build_dir, "build-$tag");
        my $windows_installer_basename = "dzlab-tools-$tag.exe";
        my $linux_installer_basename   = "dzlab-tools-$tag.linux.tar.gz";
        my $windows_installer_target   = catfile($release_dir, $windows_installer_basename);
        my $linux_installer_target     = catfile($release_dir, $linux_installer_basename);

        #######################################################################
        # windows

        if (-f $windows_installer_target && -f $linux_installer_target){
            say STDERR "$windows_installer_basename and $linux_installer_basename already exists, skipping $tag";
            next;
        }
        if (! -f $windows_installer_target){
            # clone into build dir
            in_dir($build_dir, sub {
                    warn "building windows $tag in " . getcwd();
                    remove_tree($tag_build_dir);
                    `git clone -b $tag $dzlab_tools_master $tag_build_dir `;
                });

            # checkout appropriate dir
            in_dir($tag_build_dir, sub {
                    warn "building windows $tag in " . getcwd();
                    `git checkout $tag `;
                    `tpage --define version=$tag $installer_nsi_template > installer.nsi`;
                    `tpage --define version=$tag $dzlab_check_template > dzlab-check.pl`;
                    `makensis installer.nsi`;
                    copy($windows_installer_basename, $windows_installer_target);
                });
        }
        else{
            say STDERR "$windows_installer_basename already exists";
        }

        #######################################################################
        # linux

        if (! -f $linux_installer_target){
            say STDERR "building linux $tag in " . getcwd();
            `git archive --format=tar.gz --prefix=dzlab-tools-$tag/ $tag > $linux_installer_target`;
        }
        # git archive --format=tar.gz --prefix=git-1.4.0/ v1.4.0 >git-1.4.0.tar.gz
        else{
            say STDERR "$linux_installer_target already exists";
        }

    }
}

sub tag_ok{
    my $tag_string = shift;
    my ($major, $minor, undef) = $tag_string =~ /(\d+)\.(\d+)\.(\d+)/xm;
    return looks_like_number $major 
        && looks_like_number $minor 
        && $major == 1 
        && $minor >= 5;
}

sub in_dir {
    my ($dir, $code) = @_;
    my $cwd = getcwd();
    chdir $dir;
    $code->();
    chdir $cwd;
}
