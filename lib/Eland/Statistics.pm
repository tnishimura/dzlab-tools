package Eland::Statistics;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Eland::Parser;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(eland_single_stat);
our @EXPORT = qw();

sub eland_single_stat{
    my $file_or_fh = shift;
    my %results;
    my $total_count = 0;
    my $nm_count = 0;

    my $p = Eland::Parser->new(file => $file_or_fh);
    my $counter = 0;
    while (defined(my $eland = $p->next())){
        my (undef, undef, @matches) = @$eland;
        say STDERR $counter if ++$counter % 50000 == 0;

        ++$total_count;
        ++$nm_count if ! @matches;
        for my $m (@matches) {
            my ($seq, $mismatch) = @$m;
            # there's got to be a better way of handling cases...
            ++$results{lc $seq}{$mismatch};
            ++$results{lc $seq}{'total'};
        }
    }
    return ($total_count, $nm_count, \%results);
}

# sub eland_pair_stat{
#     my ($left_file, $right_file) = @_;
# 
#     my %results;
#     open my $left_fh, '<', $left_file;
#     open my $right_fh, '<', $right_file;
# 
#     while (defined(my $right_line = <$right_fh>) && defined(my $left_line = <$left_fh>)){
#         my (undef, undef, @left_matches) = parse_eland($left_line);
#         my (undef, undef, @right_matches) = parse_eland($right_line);
#         for my $lmatch (@left_matches) {
#             my ($lseq, $lmismatch) = @$lmatch;
#             for my $rmatch (@right_matches) {
#                 my ($rseq, $rmismatch) = @$rmatch;
#                 ++$results{$lseq}{$rseq}{$lmismatch}{$rmismatch};
#             }
#         }
#     }
# 
#     close $left_fh;
#     close $right_fh;
# 
#     return \%results;
# }

1;
