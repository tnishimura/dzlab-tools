package FastqReader::CheckPairedIntegrity;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(check_paired_integrity);

my $id_line_regex = qr{^(@.*/)(\d)};

sub check_paired_integrity{
    my ($left_file_or_handle, $right_file_or_handle) = @_;
    my $left_fqr = FastqReader->new(file => $left_file_or_handle);
    my $right_fqr = FastqReader->new(file => $right_file_or_handle);

    my $left_line = 1;
    my $right_line = 1;
    while (!$left_fqr->eof() || !$right_fqr->eof()){
        my $left = $left_fqr->next();
        my $right = $right_fqr->next();

        croak "uneven lines in /1 file" if ! $left;
        croak "uneven lines in /2 file" if ! $right;

        if ($left->[0] =~ $id_line_regex){
            my $left_body = $1;
            my $left_number = $2;
            if ($right->[0] =~ $id_line_regex){
                my $right_body = $1;
                my $right_number = $2;

                if ($left_body ne $right_body || $left_number != 1 || $right_number != 2){
                    say "mismatch at $left_line in /1, $right_line in /2:";
                    say "@{$left}@{$right}";
                    exit 1;
                }
            }
        }
        $left_line += 4;
        $right_line += 4;
    }
}

1;

# 
