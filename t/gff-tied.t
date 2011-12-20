#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Test::More qw(no_plan);
use FindBin;
use lib "$FindBin::Bin/../lib";
use GFF;
use GFF::Tied;

my $tied = GFF::Tied->new(file => 't/data/gff-tied.gff');

### first 8 lines should be blank

ok(! defined $tied->get(-1), 'out of bounds before');
ok(! defined $tied->get(-999), 'out of bounds before');
ok(! defined $tied->get(18), 'out of bounds after');
ok(! defined $tied->get(999), 'out of bounds after');

for (0 .. 7){ 
    is($tied->get($_),0, "comment/blank should be zero");
}

### no attirbutes. 

{
    ok(!$tied->get(8)->equals(
            sequence => 'ctg123-meow', 
            source => undef, 
            feature => 'mRNA', 
            start => 1050, 
            end => 9000, 
            frame => undef,
            strand => '+',
            score => undef,
            attribute_string => "ID=mRNA00001;Parent=gene00001;Name=EDEN.1"
        ),
        "gff line 1- should NOT be equal");
}

{
    ok($tied->get(9)->equals(
            sequence => 'ctg123', 
            source => undef, 
            feature => 'mRNA', 
            start => 1050, 
            end => 9000, 
            frame => undef,
            strand => '+',
            score => undef,
            attribute_string => "ID=mRNA00001;Parent=gene00001;Name=EDEN.1,123",
        )
        , "gff line 2- no attr");
}

{
    ok($tied->get(10)->equals(
            sequence => 'ctg123', 
            source => undef, 
            feature => 'mRNA', 
            start => 1050, 
            end => 9000, 
            frame => undef,
            strand => '+',
            score => undef,
            attribute_string => "mRNA00001",
        ), "gff line 3- no attr");
}

### With attributes. don't pass attribute_string b/c don't care if we have parsed attrs

for (1..3){
    ok($tied->get(11)->equals(
            sequence => 'ctg123', 
            source => undef, 
            feature => 'mRNA', 
            start => 1050, 
            end => 9000, 
            frame => undef,
            strand => '+',
            score => undef,
            ID => "mRNA00001",
            Parent => "gene00001",
            Name => "EDEN.1",
        ), "gff with  attributes line 1 (try $_)");
}

{
    ok($tied->get(12)->equals(
            sequence => 'ctg123', 
            source => undef, 
            feature => 'mRNA', 
            start => 1050, 
            end => 9000, 
            frame => undef,
            strand => '+',
            score => undef,
            ID => "mRNA00001",
            Parent => "gene00001",
            Name => "EDEN.1,123",
        ), "gff with  attributes line 2");
}

{
    ok($tied->get(13)->equals(
            sequence => 'ctg123', 
            source => undef, 
            feature => 'mRNA', 
            start => 1050, 
            end => 9000, 
            frame => undef,
            strand => '+',
            score => undef,
            Note => "mRNA00001",
        ), "gff with  attributes line 3");
}

for (1..3){
    ok(!$tied->get(14)->equals(
            sequence => 'hello', 
            source => undef, 
            feature => 'exon', 
            start => 199, 
            end => 1233, 
            frame => 2,
            strand => '-',
            score => 1.2,
            c => 123,
            n => 666,
            t => 4,
        ), "gff with  attributes line 4- should NOT be equal (try $_)");
}

{
    ok($tied->get(15)->equals(
            sequence => 'hello', 
            source => undef, 
            feature => 'exon', 
            start => 199, 
            end => 1233, 
            frame => 2,
            strand => '-',
            score => 1.2,
            c => 123,
            n => 1234,
            t => 4,
        ), "gff with  attributes line 5- trailing whitespace");
}

{
    ok($tied->get(16)->equals(
            sequence => 'hello', 
            source => undef, 
            feature => 'exon', 
            start => 199, 
            end => 1233, 
            frame => 2,
            strand => '-',
            score => 1.2,
            c => "123",
            n => 1234,
            t => "4",
        ), "gff with  attributes line 5- trailing && embedded whitespace");
}

for (1..3){
    is_deeply([$tied->get(17)->parse_locus('ID')],['mRNA0001','3'], "parsing locus (try $_)");
}

