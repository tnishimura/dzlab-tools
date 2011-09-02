#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util qw /sum/;
use feature 'say';
use FindBin;
use lib "$FindBin::Bin/lib";
use DZUtil qw/open_maybe_compressed fastq_convert_read_header/;

# Check required command line parameters
pod2usage( -verbose => 1 )
    unless @ARGV;

my $type = 'verbose';
my $frequencies;
my $paired;
my $gff;
my $id_regex;
my $reference;
my $output;
my $unmatched;
my @splice;
my $no_normalize;
my $whole;

# Grabs and parses command line options
my $result = GetOptions(
    'output|o=s'    => \$output,
    'recover|u=s'   => \$unmatched,
    'splice|s=i{2}' => \@splice,
    'type|t=s'      => \$type,
    'frequencies|f' => \$frequencies,
    'paired|p'      => \$paired,
    'gff|g'         => \$gff,
    'id-regex|i=s'  => \$id_regex,
    'reference|r=s' => \$reference,
    'no-normalize|n'=> \$no_normalize,
    'whole|w'       => \$whole,
    'verbose|v'     => sub { use diagnostics; },
    'quiet|q'       => sub { no warnings; },
    'help|h'        => sub { pod2usage( -verbose => 1 ); },
    'manual|m'      => sub { pod2usage( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage( -verbose => 1 )
    unless $result
        and ( $type eq 'concise' xor $type eq 'verbose' )
        and (  ( $frequencies xor $paired )
            or ( !$frequencies and !$paired ) );

# redirect standard output to file if requested
if ($output && $output eq '-') {
    # stdout
}
elsif ($output) {
    open my $USER_OUT, '>', $output
        or croak "Can't open $output for writing: $!";
    select $USER_OUT;
} else {
    use File::Spec::Functions;
    use File::Basename;
    my ($basename, $path, undef) = fileparse($ARGV[0], ".bowtie");
    my $outfile = catfile($path,$basename) . ".eland";
    open my $USER_OUT, '>', $outfile
        or croak "Can't open $outfile for writing: $!";
    select $USER_OUT;
}

# read in bowtie verbose file
my $counts   = undef;
my $previous = undef;

while (<>) {
    chomp;
    s/[\r\n]//;

    my $current = read_bowtie($_);

    $current->{snps}{ $current->{snp}->[0] }++;

    if ($frequencies) {
        $counts->{ $current->{target}->[0] }{alternatives}
            += $current->{alternatives};
        $counts->{ $current->{target}->[0] }{frequencies}++;
    }
    elsif ($gff) {
        print_gff( $current );
    }
    elsif ($paired) {

        my $next = <>;
        chomp $next;
        $next =~ s/[\r\n]//;
        $next = read_bowtie($next);

        $next->{snps}{ $next->{snp}->[0] }++;

        print_paired_gff( $current, $next );
    }
    else {
        if ( ! defined $previous ) {
            $previous = $current;
        }
        elsif ( $current->{read_id} eq $previous->{read_id} ) {
            push @{ $previous->{strand} },     $current->{strand}->[0];
            push @{ $previous->{target} },     $current->{target}->[0];
            push @{ $previous->{coordinate} }, $current->{coordinate}->[0];
            push @{ $previous->{snp} },        $current->{snp}->[0];
            $previous->{snps}{ $current->{snp}->[0] }++;
        }
        else {

            catch_up( $previous, $unmatched, @splice )
                if $unmatched;

            print_eland($previous, $type);
            $previous = $current;
        }
    }
}

count_reads( $reference, $counts, $id_regex ) if $frequencies;

catch_up( $previous, $unmatched, @splice )
    if (defined $previous and $unmatched);

print_eland($previous, $type) if defined $previous;

# for when last read in bowtie file is *not* last read in fasta file
catch_up( $previous, $unmatched, @splice )
    if (defined $previous and $unmatched);

# T:
# Catchup: when unmatched is given, read from the original fasta file up to current read and output
# with a NM.  This is why Fasta needs to be in same order as bowtie file, and why bowtie needs to be run 
# WITHOUT multithreading
{
    my %file_handles;
    my $file_handle;

    sub catch_up {
        my ( $current, $unmatched, @splice ) = @_;


        if (! exists $file_handles{$unmatched}){
            $file_handles{$unmatched} = open_maybe_compressed($unmatched);
            #open $file_handle, '<', $unmatched
            #or croak "Can't open $unmatched: $!";
            #$file_handles{$unmatched} = $file_handle;
        }

        $file_handle = $file_handles{$unmatched};

    FASTA_HEADER:
        while ( defined( my $header = <$file_handle> )
            and defined( my $sequence = <$file_handle> ) )
        {

            chomp $header;
            chomp $sequence;

            $header =~ s/^([>@])//;
            if ( q{@} eq $1 ) {
                $header = fastq_convert_read_header($header);
                <$file_handle>;
                <$file_handle>;
            }
            elsif ( q{>} ne $1 ) {
                croak
                    "Can't figure out whether this file is straight fasta or fastq";
            }

            my %unmatched = ( header => $header, sequence => $sequence );

            # if potentially unmatched read is not current bowtie read
            #if ( $unmatched{header} !~ m/$current->{read_id}/ ) { # BUG. why regexp?
            if ( $unmatched{header} ne $current->{read_id} ) {

                $unmatched{sequence} = substr $unmatched{sequence},
                    ( $splice[0] - 1 ), ( $splice[1] - $splice[0] + 1 )
                    if (@splice && ! $whole);

                print join( "\t",
                    $unmatched{header}, $unmatched{sequence}, 'NM', "\n" );
            }
            else {
                # caught up. fetch original, splice if necessary, or fix coordinate if not.
                $current->{sequence} = $sequence;
                #say STDERR Dumper $current;

                if (@splice && $whole){
                    $current->{coordinate} = [map {
                        $_ - $splice[0] + 1;
                    } @{$current->{coordinate}}];
                } elsif (@splice && ! $whole){
                    $current->{sequence} = substr $current->{sequence},
                    ( $splice[0] - 1 ), ( $splice[1] - $splice[0] + 1 );
                }

                last FASTA_HEADER;
            }
        }
    }
}


sub print_gff {
    my ($eland) = @_;

    print join( "\t",
                $eland->{target}->[0],
                'bowtie',
                'read',
                $eland->{coordinate}->[0],
                $eland->{coordinate}->[0] + length $eland->{sequence},
                q{.},
                $eland->{strand}->[0],
                q{.},
                "read=$eland->{read_id}; alt=$eland->{alternatives}; snp=$eland->{snp}->[0]", ), "\n";
}  

sub print_paired_gff {
    my ( $current, $next ) = @_;

    print join( "\t",
        $current->{target}->[0],
        'bowtie',
        'frag',
        $current->{coordinate}->[0],
        $next->{coordinate}->[0] + length( $next->{sequence} ) - 1,
        q{.},
        q{+},
        q{.},
        "alt=$current->{alternatives}" ), "\n";
}

# eland format:
# HS2_0069:5:1:1232:2189#0/1	CGATTTTAAAAGTTTAAGTAGTGTTTTTTTGTTAGAAGATATAAAGTTAAAGATTTATATGGATTTTGGATNNATTATNNNNNNTTTGAGAAGTAGGAAG	NM	
# HS2_0069:5:1:1244:2225#0/1	TTATTTAGGATTGATAAGAATAGTTTTGAGGGAATAAATGTTAATTGATTTTTTTGTTGTTGATATTAATGNNAGTTTNNNNNNTGTGGTTTTGAATAAG	NM	
# HS2_0069:5:1:1394:2132#0/1	NCAACCGTATTTTAAAGGCGTAAGAATTGTATTCTAGTTAAAAGATATAAAGTTAAAGATTTATATGGATTTTGGCTATATTATGAAAGTTTTGAGAAGC	1:0:0:0	CHR1:15087567F0
# HS2_0069:5:1:1292:2162#0/1	NTTAATAAAAATCAGGATATATTATAAATGAAGCAGTGGTGGTTATTTTTAAGAAGGGATTGTTGTAAGGANTTTAATNNNNNNAATATTTGGCGATTAG	NM	
# HS2_0069:5:1:1353:2159#0/1	NTTTGATGGTTTAGCCGAAGTTTATGTGAGTTTTTATCTTTGTATCTTCTAACAAGGAAATATTATTTAGGTTTTTAAGATTGGTTGTGGTTTAAGTTTT	0:1:0:0	RC_CHR5:15014625F1
# HS2_0069:5:1:1265:2172#0/1	GATGAGATTTTGTGATGAAAGTATGTTGGATATAGTTGTGGAATTTAGGTATGGAATTTAGAATGTTATGTNNGTTTANNNNNNATTATATGGTGGTTAG	NM	
# HS2_0069:5:1:1373:2171#0/1	ATATATGTAACGGTTTTCGGGTTTGTTTTATTTGCAGAACTATATGAGACTGATTTTTTTATATTATATTGTTTTTGTAAGTTTTTTAAATTTTTTTGGT	1:0:0:0	CHR5:189832F0

sub print_eland {
    my ($previous, $type) = @_;

    my ( $read_id, $sequence, $chromosomes_ref, $coordinates_ref,
        $strands_ref, $mismatches_ref, $mismatches_total )
        = (
        $previous->{read_id}, $previous->{sequence},
        $previous->{target},  $previous->{coordinate},
        $previous->{strand},  $previous->{snp},
        $previous->{snps}
        );

    croak
        "Total number of chromosomes, coordinates, strands, and mismatches don't match"
        unless scalar @{$chromosomes_ref} == scalar @{$coordinates_ref}
            and scalar @{$chromosomes_ref} == scalar @{$strands_ref}
            and scalar @{$chromosomes_ref} == scalar @{$mismatches_ref};

    my $mismatches = join q{:},
        map { $mismatches_total->{$_} || 0 } ( 0 .. 3 );

    my $target;
    map {
        $target
            .= $chromosomes_ref->[$_] . q{:}
            . $coordinates_ref->[$_]
            . ( $strands_ref->[$_] eq q{+} ? q{F} : q{R} )
            . $mismatches_ref->[$_]
            . ( $_ < @{$chromosomes_ref} - 1 ? q{,} : q{} )
    } ( 0 .. @{$chromosomes_ref} - 1 );

    if ('concise' eq $type) {
        print join( "\t", $read_id, grep { $_ != 0 } split /:/, $mismatches ), "\n";
    }
    else {
        print join( "\t", $read_id, $sequence, $mismatches, $target, ), "\n";
    }
}

# HS2_0069:5:1:1742:2136#0/1      +       RC_CHR2 16657675        NAGTGGTTTGGAAGAGATGGGAAGATGGTGGGATTTGTTTTTGTTTTATA      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      0       0:A>N
# HS2_0069:5:1:1527:2151#0/1      +       RC_CHR1 17711436        NTATAGATTTTTTTTTAAATGGTGTTTATAAGTTGGTAATTTTATGGTTA      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      5       0:T>N
# HS2_0069:5:1:1606:2148#0/1      +       CHR2    13865040        NATTTTTTTTTTTTTTTTAATTTTTTTTTGTATTTTGGATTTTTTTTATT      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      0       0:T>N,8:G>T
# HS2_0069:5:1:1644:2143#0/1      +       CHRC    26379   NTTTTATTTTTTATTGAATTATATTTATAGAATTTGATTTAGTAATAATG      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      0       0:T>N


sub read_bowtie {
    my ($bowtie_line) = @_;

    my ( $read_id, $strand, $target, $coordinate, $sequence, $qualities,
        $alternatives, $snp )
    = split /\t/, $bowtie_line;
    #= split /\t|\s+/, $bowtie_line;

    # T: cause of bug! 
    #$strand = '-' if $target =~ s/^RC_//;

    # T: WTF?
    #my @mm = $snp ? split /,/, $snp : (); 
    my $mm = $snp =~ tr/>//;


    return {
        'read_id'      => $read_id,
        'strand'       => [$strand],
        'target'       => [$target],
        'coordinate'   => [$coordinate],
        'sequence'     => $sequence,
        'qualities'    => $qualities,
        'alternatives' => $alternatives,
        'snp'          => [ $mm ],
    };
}

sub count_reads {
    my ( $reference, $counts_ref, $id_regex, $no_normalize ) = @_;

    return unless $reference;

    # read in chromosome/model lengths
    my %reference = %{ index_fasta($reference) };

    $id_regex //= 'ID';

    # sort and print target id, read frequency on mapped to target per kb, average alternative mappings
  TARGET:
    for my $target ( sort keys %{$counts_ref} ) {

        my ($id) = $target =~ m/$id_regex\s*=?([^;]+)?/;

        unless ( exists $reference{$target} ) {
            carp "$target doesn't exist in $reference\n";
            next TARGET;
        }

        print join(
            "\t",

            $id,
            sprintf(
                "%g",
                $no_normalize
                ? $counts_ref->{$target}{frequencies}
                : (
                    $counts_ref->{$target}{frequencies} /
                    length $reference{$target}
                ) * 1000
            ),
            sprintf(
                "%g",
                (     $counts_ref->{$target}{frequencies}
                      ? $counts_ref->{$target}{alternatives}
                      / $counts_ref->{$target}{frequencies}
                      : 0
                  )
            ),
        ),
        "\n";
    }
}

sub index_fasta {
    my $reference_file = shift;

    my %reference = ();

    return \%reference unless $reference_file;

    # reads in the reference genome file into @fastaseq
    open my $REF, '<', "$reference_file"
        or croak "Can't open $reference for reading: $!";
    my @fastaseq = <$REF>;
    close $REF;

# find and store indices for each chromosome change and corresponding descriptions
    my ( @idx, @dsc ) = ();
    for my $i ( 0 .. @fastaseq - 1 ) {
        if ( $fastaseq[$i] =~ m/^>/ ) {
            $fastaseq[$i] =~ s/>//g;
            $fastaseq[$i] =~ s/[\r\n]//g;
            push @idx, $i;
            push @dsc, $fastaseq[$i];
        }
    }

    for my $j ( 0 .. @idx - 1 ) {
        my $line;
        if ( $j == scalar @idx - 1 ) {
            $line = join( q{}, @fastaseq[ $idx[$j] + 1 .. @fastaseq - 1 ] );
        }
        else {
            $line = join( q{},
                @fastaseq[ $idx[$j] + 1 .. $idx[ $j + 1 ] - 1 ] );
        }
        $line =~ s/[\n\r]//g;
        $reference{ $dsc[$j] } = $line;
    }
    return \%reference;
}

__END__

# 


=head1 NAME

 parse_bowtie.pl - Process bowtie alignment file into simple fields

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 OPTIONS

 parse_bowtie.pl [OPTION]... [FILE]...

 -s, --splice       splice original sequences when recovering (x y)
 -f, --frequencies  output SEQID    FREQUENCY(reads/100bp)    ALTERNATIVES
 -t, --type         type of bowtie output file (verbose or concise -- only verbose supported for now)
 -i, --id-regex     perl-type regular expression to identify feature id (ie. gene) in fasta alignment header (must include capturing parenthesis)
 -r, --reference    genome/cDNA models file in fasta format (for calculating relative frequency scores, etc.)
 -u, --recover      given original alignment fasta file, recovers unmatched reads (which bowtie does not output)
 -g, --gff          convert bowtie alignment format to gff
 -p, --paired       convert bowtie's paired ends output to gff with concatenated library ends per region
 -n, --no-normalize do not normalize frequencies
 -o, --output       filename to write results to (defaults to STDOUT)
 -v, --verbose      output perl's diagnostic and warning messages
 -q, --quiet        supress perl's diagnostic and warning messages
 -h, --help         print this information
 -m, --manual       print the plain old documentation page

=head1 REVISION

 Version 0.0.2

 $Rev: 444 $:
 $Author: psilva $:
 $Date: 2010-12-01 17:28:38 -0800 (Wed, 01 Dec 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/parse_bowtie.pl $:
 $Id: parse_bowtie.pl 444 2010-12-02 01:28:38Z psilva $:

=head1 AUTHOR

 Pedro Silva <pedros@berkeley.edu/>
 Zilberman Lab <http://dzlab.pmb.berkeley.edu/>
 Plant and Microbial Biology Department
 College of Natural Resources
 University of California, Berkeley

=head1 COPYRIGHT

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
