#!/usr/bin/perl

use warnings;
use strict;
use diagnostics;

my $loriginal=$ARGV[0];
my $roriginal=$ARGV[1];
my $mutated=$ARGV[2];
my $readsize=$ARGV[3];

die("usage: <fastafile> <fastafile> <gff align file> <read size>") unless @ARGV=4;

open(my $LORIG, "<", "$loriginal") or die("Can't open $loriginal");
open(my $RORIG, "<", "$roriginal") or die("Can't open $roriginal");
open(my $MUTATED, "<", "$mutated") or die("Can't open $mutated");


while(my $mutline = <$MUTATED>) {

#    next if((split "\t", $mutline)[5] == 0);

    if($mutline =~ m/\/1:[ACGTN]{$readsize}/) {
	my $lorigline = <$LORIG>;
	if($lorigline =~ m/^>/) {
	    $lorigline = <$LORIG>;
	}
	chomp($lorigline);
	$mutline =~ s/[ACGTN]{$readsize}/$lorigline/ ;
    }

    if($mutline =~ m/\/2:[ACGTN]{$readsize}/) {
	my $rorigline = <$RORIG>;
	if($rorigline =~ m/^>/) {
	    $rorigline = <$RORIG>;
	}
	chomp($rorigline);
	$mutline =~ s/[ACGTN]{$readsize}/$rorigline/ ;    
    }

    print $mutline;
}

close($LORIG);
close($RORIG);
close($MUTATED);
