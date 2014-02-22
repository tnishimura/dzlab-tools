#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use Data::Dumper;
use autodie;
use Test::More;
use Test::Exception;

ok(0 == system("./gff-quantile-normalization.pl t/data/gff-quantile-normalization-q1 t/data/gff-quantile-normalization-q2 t/data/gff-quantile-normalization-q3"));
ok(0 == system("head t/data/gff-quantile-normalization-q1.norm.gff t/data/gff-quantile-normalization-q2.norm.gff t/data/gff-quantile-normalization-q3.norm.gff"));

say "PLEASE COMPARE MANUALLY TO\n" . <<END;
==> q1.norm.gff <==
chr1	.	.	1	1	5.66666666666667	.	.	.
chr1	.	.	2	2	2	.	.	.
chr1	.	.	3	3	3	.	.	.
chr1	.	.	4	4	4.66666666666667	.	.	.

==> q2.norm.gff <==
chr1	.	.	1	1	4.66666666666667	.	.	.
chr1	.	.	2	2	2	.	.	.
chr1	.	.	3	3	4.66666666666667	.	.	.
chr1	.	.	4	4	3	.	.	.

==> q3.norm.gff <==
chr1	.	.	1	1	2	.	.	.
chr1	.	.	2	2	3	.	.	.
chr1	.	.	3	3	4.66666666666667	.	.	.
chr1	.	.	4	4	5.66666666666667	.	.	.
END


done_testing();
