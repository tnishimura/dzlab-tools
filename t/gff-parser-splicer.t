#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;

use FindBin;
use lib "$FindBin::Bin/../lib";

use GFF::Parser::Splicer;
use Test::GFF;
use Test::More;

my @default_colnames = qw/sequence source feature start end score strand frame attribute_string/;

my $p = GFF::Parser::Splicer->new(file => \*DATA, skip => 0, columns => \@default_colnames);

### no attirbutes. 

{
    my $gff = $p->next();

    isnt($gff->[0], 'ctg123-meow');
    ok(! defined $gff->[1]);
    is($gff->[2], 'mRNA');
    is($gff->[3], 1050);
    is($gff->[4], 9000);
    ok(! defined $gff->[5]);
    is($gff->[6], '+');
    ok(! defined $gff->[7]);
    is($gff->[8], "ID=mRNA00001;Parent=gene00001;Name=EDEN.1");
}

# {
#     my $gff = $p->next();
#     gff_has($gff,{ 
#             sequence         => 'ctg123',
#             source           => undef,
#             feature          => 'mRNA',
#             start            => 1050,
#             end              => 9000,
#             frame            => undef,
#             strand           => '+',
#             score            => undef,
#             attribute_string => "ID=mRNA00001;Parent=gene00001;Name=EDEN.1,123",
#         }, "gff line 2- no attr");
# }
# 
# {
#     my $gff = $p->next();
#     gff_has($gff,{
#             sequence         => 'ctg123',
#             source           => undef,
#             feature          => 'mRNA',
#             start            => 1050,
#             end              => 9000,
#             frame            => undef,
#             strand           => '+',
#             score            => undef,
#             attribute_string => "mRNA00001",
#         }, "gff line 3- no attr");
# }
# 
# 
# ### With attributes. don't pass attribute_string b/c don't care if we have parsed attrs
# 
# {
#     my $gff = $p->next();
#     gff_has($gff,{
#             sequence => 'ctg123',
#             source   => undef,
#             feature  => 'mRNA',
#             start    => 1050,
#             end      => 9000,
#             frame    => undef,
#             strand   => '+',
#             score    => undef,
#             ID       => "mRNA00001",
#             Parent   => "gene00001",
#             Name     => "EDEN.1",
#         }, "gff with  attributes line 1");
# }
# 
# {
#     my $gff = $p->next();
#     gff_has($gff,{
#             sequence => 'ctg123',
#             source   => undef,
#             feature  => 'mRNA',
#             start    => 1050,
#             end      => 9000,
#             frame    => undef,
#             strand   => '+',
#             score    => undef,
#             ID       => "mRNA00001",
#             Parent   => "gene00001",
#             Name     => "EDEN.1,123",
#         }, "gff with  attributes line 2");
# }
# 
# {
#     my $gff = $p->next();
#     gff_has($gff,{
#             sequence => 'ctg123',
#             source   => undef,
#             feature  => 'mRNA',
#             start    => 1050,
#             end      => 9000,
#             frame    => undef,
#             strand   => '+',
#             score    => undef,
#             Note     => "mRNA00001",
#         }, "gff with  attributes line 3");
# }
# 
# {
#     my $gff = $p->next();
#     not_gff_has($gff,{
#             sequence => 'hello',
#             source   => undef,
#             feature  => 'exon',
#             start    => 199,
#             end      => 1233,
#             frame    => 2,
#             strand   => '-',
#             score    => 1.2,
#             c        => 123,
#             n        => 666,
#             t        => 4,
#         }, "gff with  attributes line 4- should NOT be equal");
# }
# 
# {
#     my $gff = $p->next();
#     gff_has($gff,{
#             sequence => 'hello',
#             source   => undef,
#             feature  => 'exon',
#             start    => 199,
#             end      => 1233,
#             frame    => 2,
#             strand   => '-',
#             score    => 1.2,
#             c        => 123,
#             n        => 1234,
#             t        => 4,
#         }, "gff with  attributes line 5- trailing whitespace");
# }
# 
# {
#     my $gff = $p->next();
#     gff_has($gff,{
#             sequence => 'hello',
#             source   => undef,
#             feature  => 'exon',
#             start    => 199,
#             end      => 1233,
#             frame    => 2,
#             strand   => '-',
#             score    => 1.2,
#             c        => "123",
#             n        => 1234,
#             t        => "4",
#         }, "gff with  attributes line 5- trailing && embedded whitespace");
# }
# 
# {
#     my $gff = $p->next();
#     is_deeply([$gff->parse_locus('ID')],['mRNA0001','3'], "parsing locus");
# }

done_testing();



__DATA__

 
#
 #
 # 
#helo
 #helo
 #helo asd 
## pragma 1
 ## pragma 2
 ## 
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1,123
ctg123	.	mRNA	1050	9000	.	+	.	mRNA00001
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1,123
ctg123	.	mRNA	1050	9000	.	+	.	mRNA00001
hello	.	exon	199	1233	1.2	-	2	c=123;n=1234;t=4
hello	.	exon	199	1233	1.2	-	2	c=123;n=1234;t=4  
hello	.	exon	199	1233	1.2	-	2	  c  =  123  ;  n  =  1234  ;  t  =  4  
hello	.	exon	199	1233	1.2	-	2	ID= mRNA0001.3 ;  c  =  123  ;  n  =  1234  ;  t  =  4  
