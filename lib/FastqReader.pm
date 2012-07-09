package FastqReader;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Moose;
use Carp;
use autodie;    
use List::MoreUtils qw/ any /;
use DZUtil qw/c2t g2a/;

with 'ParserRole';

has 'linesper' => (
    is => 'rw',
    default => 4,
);

# given a list of read ids, return a hash of { id => sequence }
sub get_reads{
    my $self = shift;
    my @query_ids = @_;
    return {} if ! @query_ids;

    my %results;
    while (defined(my $q = $self->next())){
        my ($readid, $sequence) = @$q;
        QUERY:
        for my $query (@query_ids) {
            if ($readid =~ /\Q$query\E/){
                $results{$query} = $sequence;
                last QUERY;
            }
        }
    }
    return \%results;
}

sub next{
    my $self = shift;
    my $fh = $self->filehandle;

    defined(my $readid = scalar readline $fh) or return;
    defined(my $sequence = scalar readline $fh) or return;
    chomp($readid, $sequence);

    if ($self->linesper() == 4){
        defined(my $readid_again = scalar readline $fh) or return;
        defined(my $quality = scalar readline $fh) or return;
        chomp($readid_again, $quality);

        $readid =~ s/\^@//;
        $readid_again =~ s/\^\+//;
        return [$readid, $sequence, $readid_again, $quality];
    }
    elsif ($self->linesper()==2){
        $readid =~ s/\^>//;
        return [$readid, $sequence];
    }
    else{
        croak "only linesper 2,4 supported";
    }
}

sub fastq_to_fasta{
    my ($methyl_pattern, $input_file_or_filehandle, $output_file_or_filehandle) = @_;
    my $fh;
    if (ref $output_file_or_filehandle eq 'GLOB'){
        $fh = $output_file_or_filehandle;
    }
    else{
        open my $tmpfh, '>', $output_file_or_filehandle;
        $fh = $tmpfh;
    }
    if (defined($methyl_pattern) and ($methyl_pattern ne 'c2t' and $methyl_pattern ne 'g2a')){
        croak "bs_fastq expects methyl_pattern to be c2t or g2a";
    }

    my $fqr = FastqReader->new(file => $input_file_or_filehandle);
    while (defined(my $quartet = $fqr->next())){
        say $fh ">$quartet->[0]";
        if (defined $methyl_pattern){
            if ($methyl_pattern eq 'c2t'){
                say $fh c2t($quartet->[1]);
            }
            else{
                say $fh g2a($quartet->[1]);
            }
        }
        else{
            say $fh $quartet->[1];
        }
    }

    if (ref $output_file_or_filehandle ne 'GLOB'){
        close $fh;
    }
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
