#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;
use List::MoreUtils qw/all/;
use Scalar::Util qw/refaddr/;
use Pod::Usage;
use Tie::IxHash;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;

pod2usage(-verbose => 2, -noperldoc => 1) if (!@ARGV);  

tie my %files, 'Tie::IxHash'; 
%files = map { $_ => 0 } @ARGV;

# make sure line_count is the same in every file
my $line_count = do { 
    my @f = keys %files;
    my $first_file_lc = line_count(shift @f);
    if (! all { line_count($_) } @f){
        die "all files should have the same line count";
    }
    $first_file_lc;
};

# get ranked lines for each file
for my $f (keys %files) {
    $files{$f} = rank_scores($f);
}
# dumpf(\%files);

my %output_files = calculate_rank_quantiles($line_count, \%files);
# dumpf(\%files);

# blit
for my $f (keys %output_files) {
    my @gffs = @{$output_files{$f}};
    $f =~ s/\.gff$//;
    $f .= ".norm.gff";
    open my $fh, '>', $f;
    say $fh $_ for @gffs;
    close $fh;
}

#######################################################################

# accept { FILE => [GFF, SCORE, LINE, RANK] } 
# return { FILE => [GFF] }
sub calculate_rank_quantiles{
    my $line_count = shift;
    my $files = shift;
    my $file_count = scalar keys %$files;
    my @sums = (0) x $line_count;

    for my $f (keys %$files) {
        my @ranked_scores = sort { $a<=>$b} map { $_->[1] } @{$files{$f}};
        for my $i (0 .. $line_count - 1) {
            $sums[$i] += $ranked_scores[$i];
        }
    }

    # say Dumper $line_count, \@sums;
    my @ranked_quantiles = map { $_ / $file_count } @sums;
    # say Dumper \@ranked_quantiles;

    my %rv;
    for my $f (keys %$files) {
        $rv{$f} = [];
        for my $line (@{$files{$f}}) {
            my ($gff, undef, undef, $rank) = @{$line};
            $gff->score($ranked_quantiles[$rank]);
            push @{$rv{$f}}, $gff;
            # say "+ " . refaddr $gff;
        }
    }
    return %rv;
}

# returns [GFF_OBJECT, SCORE, LINE_NUM, RANK] in original order
sub rank_scores{
    my $file_name = shift;
    my $parser = GFF::Parser->new(file => $file_name, normalize => 0);
    my $line = 0;
    my @ranked_lines; # [GFF, SCORE, LINE, RANK]
    while (defined(my $gff = $parser->next())){
        push @ranked_lines, [$gff, $gff->score, $line];
        $line++;
    }
    @ranked_lines = sort { $a->[1] <=> $b->[1] } @ranked_lines;

    # if there is a tie, give them all the high rank possible
    my $last_idx;
    my $last_val;
    for my $i (0 .. $#ranked_lines) {
        my $val = $ranked_lines[$i][1];
        if (defined $last_val && $last_val == $val){
            $ranked_lines[$i][3] = $last_idx;
        }
        else{
            $ranked_lines[$i][3] = $i;
            $last_idx = $i;
            $last_val = $val;
        }
    }

    @ranked_lines = sort { $a->[2] <=> $b->[2] } @ranked_lines;

    return \@ranked_lines;
}

sub line_count {
    my $file_name = shift;
    open my $fh, '<:crlf', $file_name;
    my $count = 0;
    while (defined(my $line = <$fh>)){
        $count++;
    }
    close $fh;
    return $count;
}

sub dumpf{
    my $files = shift;
    for my $f (keys %$files) {
        $files->{$f} = rank_scores($f);
        say $f;
        for my $l (@{$files->{$f}}){
            say refaddr($l->[0]) . "\t$l->[2]: $l->[0]\trank $l->[3]";
        }
    }
}

=head1 gff-quantile-normalization.pl 

 gff-quantile-normalization.pl file1.gff file2.gff file3.gff ....

Creates file1.norm.gff file2.norm.gff file3.norm.gff with quantile-normalized
scores: http://en.wikipedia.org/wiki/Quantile_normalization

All gff files must have the same number of lines.  This script only looks at
column 6 (scores), so it probably does not make sense to have multiple
chromosome/sequences in a file.

=cut

__DATA__
From http://en.wikipedia.org/wiki/Quantile_normalization
A quick illustration of such normalizing on a very small dataset:

Arrays 1 to 3, genes A to D

A    5    4    3
B    2    1    4
C    3    4    6
D    4    2    8
For each column determine a rank from lowest to highest and assign number i-iv

A    iv    iii   i
B    i     i     ii
C    ii    iii   iii
D    iii   ii    iv
These rank values are set aside to use later. Go back to the first set of data. Rearrange that first set of column values so each column is in order going lowest to highest value. (First column consists of 5,2,3,4. This is rearranged to 2,3,4,5. Second Column 4,1,4,2 is rearranged to 1,2,4,4, and column 3 consisting of 3,4,6,8 stays the same because it is already in order from lowest to highest value.) The result is:

A    5    4    3    becomes A 2 1 3
B    2    1    4    becomes B 3 2 4
C    3    4    6    becomes C 4 4 6
D    4    2    8    becomes D 5 4 8
Now find the mean for each row to determine the ranks

A (2 1 3)/3 = 2.00 = rank i
B (3 2 4)/3 = 3.00 = rank ii
C (4 4 6)/3 = 4.67 = rank iii
D (5 4 8)/3 = 5.67 = rank iv
Now take the ranking order and substitute in new values

A    iv    iii   i
B    i     i     ii
C    ii    iii   iii
D    iii   ii    iv
becomes:

A    5.67    4.67    2.00
B    2.00    2.00    3.00
C    3.00    4.67    4.67
D    4.67    3.00    5.67
This is the new normalized values. The new values have the same distribution and can now be easily compared.
