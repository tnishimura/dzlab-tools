#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use Data::Dumper;
use autodie;
use Test::More;
use Test::Exception;
use GFF::FileInMemory;

my $gfile = GFF::FileInMemory->new(file => "t/data/test1.gff");

is($gfile->count, 99);
is_deeply([$gfile->sequences_in_file], ['Chr1', 'Chr2']);

is(-1, GFF::FileInMemory::_sorter_double('Chr1', 121124, 'Chr2', 11111));
is(-1, GFF::FileInMemory::_sorter_double('Chr1', 121124, 'Chr1', 121125));
is(1, GFF::FileInMemory::_sorter_double('Chr1', 121125, 'Chr1', 121124));
is(0, GFF::FileInMemory::_sorter_double('Chr1', 121125, 'Chr1', 121125));

my $found_gff = $gfile->find_by_seq_start('Chr1', 121124);
is($found_gff->sequence, 'Chr1');
is($found_gff->start, 121124);


# Chr1	TAIR8	gene	11649	13714	.	-	.	ID=AT1G01030.1;Name=AT1G01030;Note=NGA3 (NGATHA3),transcription factor
done_testing();

