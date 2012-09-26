#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More;
use Test::Exception;
use ParserNoMoose;
use IO::Scalar;

# handle
{
    my @lines    = qw{asdhka sdlkahl2kgl6 k234g ql2k34gt6 ql2k345tg ql};
    my $contents = join "\n", @lines;

    my $parser = ParserNoMoose->new(file => \$contents);

    while (defined(my $line = readline($parser->{handle}))){
        chomp $line;
        is($line, shift(@lines));
    }
}

done_testing;
