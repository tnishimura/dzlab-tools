#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Benchmark qw/timethese cmpthese/;
use TestUtils qw/test_gff/;
use GFF::Parser;
use GFF::Parser::XS;
use GFF::Parser::Simple;
use GFF::Parser::Splicer;
use GFF::Parser::ClosureSplicer;

# negative means seconds to run for each
# positive means number of times to run each
# 0 means run 3 seconds each

say "creating test_gff";
my $test_gff = test_gff(100_000);
say $test_gff;

say "running:";
cmpthese(5,{ 
        native => sub{
            say "native";
            my $p = GFF::Parser->new(file => $test_gff);
            my $count = 0;
            while (defined(my $gff = $p->next())){
                $count += $gff->get_column('n');
                $count += $gff->get_column('c');
                $count += $gff->get_column('t');
            }
        },
        splicer => sub{
            say "splicer";
            my $p = GFF::Parser::Splicer->new(file => $test_gff, columns => [qw/seq start end c t n/]);
            my $count = 0;
            while (1){
                my @gff = $p->next();
                #say Dumper \@gff;
                last if (0 == @gff);
                my ($seq, $start, $end, $c, $t, $n) = @gff;
                $count += $c;
                $count += $t;
                $count += $n;
            }
        },
    });


