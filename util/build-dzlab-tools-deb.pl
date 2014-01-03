#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;

if (0 != system("which fpm > /dev/null")){
    say "can't find fpm in \$PATH?";
    exit 1;
}

my @files = grep { /^(bin|conf|lib|share)/ } split /\n/, `git ls-files`;

my @git_tree_changes = grep { ! /\?\?/ && /\s(bin|conf|lib|share)/} split /\n/, `git status --porcelain`;
if (@git_tree_changes){
    say "*** git tree does not look clean?";
    say "(\n@git_tree_changes\n)" ;
    exit 1;
}

system("fpm -f -s dir -t deb -n dzlab-tools -v 1.5.56 --prefix /opt/dzlab-tools @files");
