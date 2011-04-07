package Fasta;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Carp;
use Data::Dumper;
use feature 'say';
use autodie;


require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(convert_file convert bisulfite_convert);
our @EXPORT = qw(fasta_subseq slurp_fasta format_fasta count_fasta count_fasta_complete dump_fasta);

=head1 EXPORTED FUNCTIONS

=head2 slurp_fasta('text.fasta', normalize => 1);

return hash of seqid (the first word after the '>') to the sequence. 
All sequence names upper cased if normalize => 1 (default)

=cut
sub slurp_fasta {
    my ($file, %opt) = @_;
    return {} unless $file;

    $opt{normalize} //= 1;

    my %accum = ();

    open my $fh, '<', $file or croak "Can't open $file: $?";

    my $current;
    my @buffer;
    while (defined(my $line = <$fh>)) {
        $line =~ tr/\r\n//d;
        if ($line =~ /^>(\w+)/){
            if ($current){
                $accum{$current} = join q{}, @buffer;
                @buffer = ();
            }
            $current = $opt{normalize} ? uc $1 : $1;
            
        } 
        else{
            chomp $line;
            push @buffer, $opt{normalize} ? uc $line : $line;
        }
    }

    $accum{$current} = join q{}, @buffer if ($current) ; # last seq
    $accum{$current} =~ s/\s*//g;

    close $fh;

    return \%accum;
}

=head2 fasta_get($fasta, $seqid, $start, $end, coord => 'f', rc => 1, base => 0, normalize => 1)

If coord => 'r', use coordinates with respect to the 3' end. 
If rc => 1 (default for reverse coord), reverse complement. 
If Normalize => 1 (default) means upper case seqid.

=cut

sub fasta_subseq{
    my ($fasta, $seqid, $start, $end, %opt) = @_;
    my $coord     = $opt{coord} // 'f';
    my $rc        = $opt{rc} // ($coord eq 'r');
    my $base      = $opt{base} // 1;
    my $normalize = $opt{normalize} // 1;
    $seqid = $normalize ? uc $seqid : $seqid;

    if ($end < $start){
        croak "fasta_get: end ($end) < start ($start)?";
    }

    # everything in base 0 coord now.
    $start -= $base;
    $end   -= $base;
    my $sublen = ($end-$start)+1;

    if(my $seq = $fasta->{$seqid}){
        my $totlen = length $seq;
        my $lastindex = $totlen - 1;
        my $sub;
        if ($coord eq 'f'){
            my $left = $start;
            my $right = $left + $sublen - 1; 

            croak "($start,$end) (left = $left, right = $right) out of bounds" if ($left < 0 || $right > $lastindex);
            $sub = substr($seq,$left, $sublen);
        } 
        elsif ($coord eq 'r') {
            my $left = $totlen - $end - 1; 
            my $right = $left + $sublen - 1;
            croak "($start,$end) (left = $left, right = $right) on reverse out of bounds" if ($left < 0 || $right > $lastindex);

            $sub = substr($seq,$left, $sublen);
        }
        else {
            croak "coord needs to be 'f' or 'r'";
        }
        if ($rc){
            $sub =~ tr/acgtACGT/tgcaTGCA/;
            $sub = reverse $sub;
        }
        return $sub;
    } else {
        croak "no such seqid $seqid";
    }
}


=head2 convert($seqs, '[c2t|g2a|rc]')

=cut

# convert c2t/g2a/rc in-place
sub convert{
    my ($seqs, $pattern) = @_;
    for my $s (keys %$seqs){
        $seqs->{$s} =~ tr/cC/tT/ if $pattern eq 'c2t';
        $seqs->{$s} =~ tr/gG/aA/ if $pattern eq 'g2a';
        $seqs->{$s} =~ tr/acgtACGT/tgcaTGCA/ if $pattern eq 'rc';
    }
}

=head2 convert_file('in.fasta', 'out.fasta', '[c2t|g2a|rc]'

=cut

sub convert_file{
    my ($infile,$outfile,$pattern) = @_;
    open my $in, '<', $infile;
    open my $out, '>', $outfile;

    while(my $line = <$in>) {
        if($line =~ m/^[ACGTN]+/i) {
            $line =~ tr/Cc/Tt/ if $pattern eq 'c2t';
            $line =~ tr/Gg/Aa/ if $pattern eq 'g2a';
            if ($pattern eq 'rc'){
                $line =~ tr/acgtACGT/tgcaTGCA/;
                $line = reverse $line;
            }
        }
        print $out $line;
    }
    close($in);
    close($out);
}

# bisulfite convert
# c2t of forward strand + c2t of original reverse complemented
sub bisulfite_convert{
    my ($file,$outfile) = @_;
    my $forward = slurp_fasta($file,normalize => 1);
    my $reverse;
    for my $s (keys %$forward){
        $reverse->{"RC_$s"} = reverse $forward->{$s};
        $reverse->{"RC_$s"} =~ tr/acgtACGT/tgcaTGCA/;
    }
    #say Dumper $forward;
    #say Dumper $reverse;

    convert($forward,'c2t');
    convert($reverse,'c2t');
    my @seqs = sort keys %$forward;

    my $fh;
    if ($outfile){
        open $fh, '>', $outfile;
    } else {
        $fh = \*STDOUT;
    }

    for my $s (sort @seqs){
        print $fh format_fasta($s,$forward->{$s});
        print $fh format_fasta("RC_$s",$reverse->{"RC_$s"});
    }

    if ($outfile){
        close $fh;
    }
}

sub fastq2fasta_file {
    my ($infile,$outfile) = @_;
    open my $i, '>', $infile;
    open my $o, '>', $outfile;
    
    while (<$i>) {
        if (/^@(\S+)/) {
            print $o ">$1\n";
            $_ = <>; print $i $_;
            <>; <>;
        }
    }
    close $i;
    close $o;
}

sub count_fasta_complete {
    my ($file) = @_;
    return {} unless $file;

    my %accum = ();

    open my $fh, '<', $file or croak "Can't open $file: $?";

    my $current;
    while (defined(my $line = <$fh>)) {
        $line =~ tr/\r\n//d;
        if ($line =~ /^>(\w+)/){
            $current = lc $1;
        } 
        else{
            my ($len, $bp) = (length $line, $line =~ tr/acgtACGT//);
            $accum{$current}{length} += $len;
            $accum{$current}{bp}     += $bp;
            $accum{$current}{nonbp}  += $len-$bp;
        }
    }
    close $fh;

    return \%accum;
}

=head2 count_fasta('file')

Return hash of {seqid => count}.

=cut
sub count_fasta {
    my ($file) = @_;
    my $counts = count_fasta_complete($file);
    return {map { $_ => $counts->{$_}{length} } keys %$counts};
}


=head2 format_fasta('header', $seq)

format a header and a single seq for printing as a fasta.

=cut 
sub format_fasta{
    my ($header, $seq, $width) = @_;

    confess "need a sequence and header" unless ($seq and $header);

    my @buffer;
    $buffer[0] = ">$header" if defined $header;
    $width //= 80;
    my $length = length $seq;

    for (my $position = 0; $position < $length; $position += $width){
        push @buffer, substr $seq, $position, $width;
    }
    return (join "\n", @buffer) . "\n";
}

=head2 dump_fasta($seqs)

Print fasta to STDOUT.

=cut

sub dump_fasta{
    my ($fasta) = @_;

    for my $seq (sort keys %$fasta) {
        say format_fasta(uc $seq, $fasta->{$seq});
    }
}

1;


