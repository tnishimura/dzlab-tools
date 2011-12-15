#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Cwd qw/getcwd/;
use File::Basename qw/basename dirname/;
use File::Path qw/make_path remove_tree/;
use File::Spec::Functions qw/rel2abs canonpath catdir catfile updir/;
use File::Copy;
use File::Temp qw/tempdir tempfile/;
use Test::More qw(no_plan);
use Perl6::Slurp;
use Launch qw/cast/;

require_ok("generic_mashup.pl");

is_deeply(
    parse_nicknames([qw/file1.txt file2.txt file3.txt/], "nick1,nick2,nick3"), 
    [["file1.txt" => 'nick1'], ["file2.txt" => 'nick2'], ['file3.txt' => 'nick3']], 
    "parse_nicknames",
);

is_deeply(
    parse_key_specs(" 1 , 2, 4, 5n"),
    [
    [1, 't'],
    [2, 't'],
    [4, 't'],
    [5, 'n'],
    ],
    "parse_key_specs 1"
);

is_deeply(
    parse_key_specs(" 2n, 4, 1, 5n"),
    [
    [2, 'n'],
    [4, 't'],
    [1, 't'],
    [5, 'n'],
    ],
    "parse_key_specs 2"
);


my $cmp_name = "cmp01";
{
    my $cmp = make_comparator( parse_key_specs("1,2n") );

    is($cmp->(['chr1' , 50]  , ['chr1' , 100]) , -1 , $cmp_name++);
    is($cmp->(['chr1' , 100] , ['chr1' , 100]) , 0  , $cmp_name++);
    is($cmp->(['chr1' , 100] , ['chr1' , 50])  , 1  , $cmp_name++);

    is($cmp->(['chr1' , 50]  , ['chr2' , 100]) , -1 , $cmp_name++);
    is($cmp->(['chr1' , 100] , ['chr2' , 100]) , -1 , $cmp_name++);
    is($cmp->(['chr1' , 100] , ['chr2' , 50])  , -1 , $cmp_name++);

    is($cmp->(['chr2' , 50]  , ['chr1' , 100]) , 1  , $cmp_name++);
    is($cmp->(['chr2' , 100] , ['chr1' , 100]) , 1  , $cmp_name++);
    is($cmp->(['chr2' , 100] , ['chr1' , 50])  , 1  , $cmp_name++);
}

{
    my $cmp = make_comparator( parse_key_specs("3n,1") );
    my $junk = "junk01";

    is($cmp->([ 50  ,'chr1' ]  , [100, 'chr1',]) , -1 , $cmp_name++);
    is($cmp->([ 100 ,'chr1' ] , [100, 'chr1',]) , 0  , $cmp_name++);
    is($cmp->([ 100 ,'chr1' ] , [50, 'chr1',])  , 1  , $cmp_name++);
                            
    is($cmp->([ 50  ,'chr1' ]  , [100, 'chr2',]) , -1 , $cmp_name++);
    is($cmp->([ 100 ,'chr1' ] , [100, 'chr2',]) , -1 , $cmp_name++);
    is($cmp->([ 100 ,'chr1' ] , [50, 'chr2',])  , 1 , $cmp_name++);
                            
    is($cmp->([ 50  ,'chr2' ]  , [100, 'chr1',]) , -1  , $cmp_name++);
    is($cmp->([ 100 ,'chr2' ] , [100, 'chr1',]) , 1  , $cmp_name++);
    is($cmp->([ 100 ,'chr2' ] , [50, 'chr1',])  , 1  , $cmp_name++);
}


{
    my $file1_sorted = "t/data/generic_mashup_sorted.txt";
    my $file1_unsorted = "t/data/generic_mashup_unsorted.txt";

    die "check data files" unless grep { -f } $file1_sorted, $file1_unsorted;

    my (undef, $tempfile) = tempfile(UNLINK => 1);
    copy($file1_unsorted, $tempfile) or die "can't copy file?";

    sort_file_in_place($tempfile, parse_key_specs("1,2n"));

    is(slurp($tempfile), slurp($file1_sorted), "sort_file_in_place");
}

is_deeply([extract_columns([qw/a b c d e/], [1,2])], [qw/a b/], "extract_columns");
is_deeply([extract_columns([qw/a b c d e/], [2,1])], [qw/b a/], "extract_columns");
is_deeply([extract_columns([qw/a b c d e/], [3,1])], [qw/c a/], "extract_columns");

{
    my $tempdir = tempdir(CLEANUP => 1);
    my $csv1 = "t/data/generic_mashup_csv1.txt";
    my $csv2 = "t/data/generic_mashup_csv2.txt";
    my $result = "t/data/generic_mashup_csv_results.txt";

    my $tmp1 = catfile($tempdir, basename($csv1));
    my $tmp2 = catfile($tempdir, basename($csv2));
    my $tmpresult = catfile($tempdir, "result");

    copy($csv1, $tmp1);
    copy($csv2, $tmp2);

    ok(cast("perl", "generic_mashup.pl", "-o", $tmpresult, '-v', '3', '-k', '1,2n', $tmp1, $tmp2), "generic_mashup.pl run");
    is(slurp($tmpresult), slurp($result), "generic_mashup.pl results ok");
}
