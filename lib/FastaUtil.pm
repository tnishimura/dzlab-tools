package FastaUtil;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use FastaReader;
use DZUtil qw/c2t g2a/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw( rc_fasta_on_disk bs_fasta_on_disk bsrc_fasta_on_disk);

=head1 USAGE

 rc_fasta_on_disk ($input_file_or_fh, $output_file_or_fh)

 bs_fasta_on_disk ('c2t' | 'g2a', $input_file_or_fh, $output_file_or_fh)

 bsrc_fasta_on_disk ('c2t' | 'g2a', $input_file_or_fh, $output_file_or_fh)

=cut
sub bs_fasta_on_disk{
    my ($pattern, $input_file_or_filehandle, $output_file_or_filehandle) = @_;
    
    if ($pattern ne 'c2t' && $pattern ne 'g2a'){
        croak "bs_convert_on_disk invalid pattern";
    }

    my $input_fh;
    my $output_fh;
    if (ref $input_file_or_filehandle eq 'GLOB'){
        $input_fh = $input_file_or_filehandle;
    }
    else{
        open $input_fh, '<', $input_file_or_filehandle;
    }

    if (ref $output_file_or_filehandle eq 'GLOB'){
        $output_fh = $output_file_or_filehandle;
    }
    else{
        open $output_fh, '>', $output_file_or_filehandle;
    }

    while(my $line = <$input_fh>) {
        if($line !~ m/^>/i) {
            $line = c2t($line) if $pattern eq 'c2t';
            $line = g2a($line) if $pattern eq 'g2a';
        }
        print $output_fh $line;
    }

    if (ref $input_file_or_filehandle ne 'GLOB'){
        close $input_fh;
    }
    if (ref $output_file_or_filehandle ne 'GLOB'){
        close $output_fh;
    }
}

sub bsrc_fasta_on_disk{
    my ($pattern, $in, $out) = @_;

    # pattern can be undef, meaning no bs-ing: for implementation of
    # rc_fasta_on_disk
    my @bsarg;
    if (! defined $pattern){
        @bsarg = ();
    }
    elsif ($pattern eq 'c2t' || $pattern eq 'g2a'){
        @bsarg = (bs => $pattern);
    }
    else{
        croak "bsrc_fasta_on_disk invalid pattern";
    }


    my $f = FastaReader->new(file => $in, slurp => 0);
    my $outfh;
    if (ref $out eq 'GLOB'){
        $outfh = $out;
    }
    else {
        open $outfh, '>', $out;
    }

    for my $seq ($f->sequence_list()) {
        #say $outfh ">$seq";
        print $outfh $f->get_pretty($seq, $seq, undef, undef, @bsarg);

        #say $outfh ">RC_$seq";
        print $outfh $f->get_pretty("RC_$seq", $seq, undef, undef, rc => 1, @bsarg);
    }
    if (ref $out ne 'GLOB'){
        close $outfh;
    }
}

sub rc_fasta_on_disk{
    my ($in, $out) = @_;
    return bsrc_fasta_on_disk(undef, $in, $out);
}

1;

