package Run::Bowtie;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use List::MoreUtils qw/any notall/;
use DZUtil qw/fastq_read_length/;
use IPC::Cmd qw/can_run run/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(bowtie);

sub _construct_common_args{
    my %opt = @_;
    my @args;
    if (! exists $opt{index} || notall { -s "$opt{index}.$_" } qw/ 1.ebwt 2.ebwt 3.ebwt 4.ebwt rev.1.ebwt rev.2.ebwt/){
        croak "need index";
    }
    push @args, $opt{index};

    # maxhits
    if (defined $opt{maxhits}  and $opt{maxhits} > 0){
        if (any {exists $opt{$_}} qw/suppress reportmax strata best/){
            croak "suppress, reportmax, strata, best are set automatically with maxhits";
        }
        push @args, -k => $opt{maxhits}, -m => $opt{maxhits}, '--strata', '--best';
    }
    else{
        push @args, (-k => $opt{reportmax}) if exists $opt{reportmax} && $opt{reportmax} >= 1;
        push @args, (-m => $opt{suppress})  if exists $opt{suppress} && $opt{suppress} >= 1;
        push @args, '--strata'              if exists $opt{strata} && $opt{strata} == 1;
        push @args, '--best'                if exists $opt{best} && $opt{best} == 1;
    }
    # splice
    if (exists $opt{splice}){
        if (any {exists $opt{$_}} qw/trim5 trim3/){
            croak "trim5, trim3 are set automatically with splice";
        }
        if (! exists $opt{readlength}){
            $opt{readlength} = fastq_read_length($opt{'-1'});
            warn "splice needs readlength option. inferring from $opt{-1}.";
        }
        my ($left, $right) = @{$opt{splice}};
        push @args, (-5 => $left - 1, -3 => $opt{readlength} - $right);
    }
    else{
        push @args, (-5 => $opt{trim5}) if exists $opt{trim5} && $opt{trim5} >= 1;
        push @args, (-3 => $opt{trim3}) if exists $opt{trim3} && $opt{trim3} >= 1;
    }
    
    push @args, (-B => $opt{base} // 1);
    push @args, (-v => $opt{mismatches} // 2);

    push @args, '--norc'                 if exists $opt{norc} && $opt{norc} == 1;
    push @args, (-p => $opt{threads})    if exists $opt{threads} && $opt{threads} == 1;
    push @args, ('--seed' => $opt{seed}) if exists $opt{seed};

    given ($opt{format}){
        when ('fasta') { push @args, '-f'; }
        when ('fastq') { push @args, '-q'; }
        when ('raw')   { push @args, '-r'; }
        when (undef)   { push @args, '-q'; }
        default        { push @args, '-q'; }
    }

    return @args;
}


=head2 bowtie

Most basic bowtie invocation. Read from file, output to file, log returned

    my ($processed, $aligned, $suppressed, $reported, @loglines) = bowtie(
        '-1'       => $reads,
        output     => $output,
        index      => $ref,
        maxhits    => 10,
        splice     => [0,50],
        readlength => 100,
    );

=cut
sub bowtie{
    croak "bowtie: argument error, uneven" if (@_ % 2 != 0);
    if (! can_run('bowtie')){
        croak "bowtie not installed?";
    }

    my %opt = @_;
    my @args = _construct_common_args(%opt);

    my $verbose = exists $opt{verbose} && $opt{verbose} >= 1;

    if (notall { exists $opt{$_} } qw/-1 output/){
        croak "need -1, output";
    }

    if (exists $opt{-1}){
        croak "-1 unreadable" if ! -s $opt{-1};
        if (exists $opt{-2}){
            push @args, ('-1', $opt{-1});
        }
        else{
            push @args, $opt{-1};
        }
    }
    if (exists $opt{-2}){
        croak "-2 unreadable" if ! -s $opt{-2};
        croak "-2 specified without -1?" if ! -s $opt{-1};
        push @args, ('-2', $opt{-2});
    }
    push @args, $opt{output};
    
    #######################################################################
    # run 
    my $bowtie_cmd = join " ", 'bowtie', @args;
    
    say STDERR $bowtie_cmd if $verbose;
    my $output_buffer;
    if (scalar run( command => $bowtie_cmd,
                    verbose => 0,
                    buffer  => \$output_buffer,
                    timeout => 20 )
    ) {
        say $output_buffer;
        my @split = split /\n/, $output_buffer;
        return _parse_bowtie_log(@split), @split;
    }
    else{
        croak "running bowtie with IPC::Cmd failed";
    }
}

sub bowtie_pipe{

    croak "bowtie: argument error, uneven" if (@_ % 2 != 0);

    my %opt = @_;
    my @args = _construct_common_args(%opt);

    my $verbose = exists $opt{verbose} && $opt{verbose} >= 1;

    if (any {exists $opt{$_}} qw/-1 -2/){
        croak "can't specify input files with bowtie_pipe";
    }
    push @args, '-';

    if (! exists $opt{output}){
        croak "need output";
    }
    push @args, $opt{output};
    
    #######################################################################
    # run 
    
    say STDERR join " ", 'bowtie', @args if $verbose;

    my $pid = open my $bowtie_process, '|-';
    defined $pid or croak "couldn't fork!";

    my @log;
    if ($pid){ # parent
        while (defined(my $logline = <$bowtie_process>)){
            chomp $logline;
            push @log, $logline;
        }
        {
            no autodie qw/close/;
            close $bowtie_process || croak "bowtie ended prematurely?";
        }
    } 
    else{
        close STDERR;
        open STDERR, ">&STDOUT";
        exec 'bowtie', @args;
    }
    return _parse_bowtie_log(@log), @log;
}


sub _parse_bowtie_log{
    my @loglines = @_;
    
    my $processed;
    my $aligned;
    my $reported;
    my $suppressed;
    if ($loglines[0] =~ /reads\sprocessed:\s(\d+)/){
        $processed = $1;
    }
    if ($loglines[1] =~ /reads\swith\sat\sleast\sone\sreported\salignment:\s(\d+)/){
        $aligned = $1;
    }
    if ($loglines[-2] =~ /reads\swith\salignments\ssuppressed\sdue\sto\s-m:\s(\d+)/){
        $suppressed = $1;
    }
    else{
        $suppressed = 0;
    }
    if ($loglines[-1] =~ /Reported\s(\d+)\salignments\sto\s\d+\soutput\sstream/){
        $reported = $1;
    }

    return ($processed, $aligned, $suppressed, $reported);
}

1;

__END__
# reads processed: 5100
# reads with at least one reported alignment: 3 (0.06%)
# reads that failed to align: 5097 (99.94%)
Reported 3 alignments to 1 output stream(s)

Usage: 
  bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]

  <m1>    Comma-separated list of files containing upstream mates (or the
          sequences themselves, if -c is set) paired with mates in <m2>
  <m2>    Comma-separated list of files containing downstream mates (or the
          sequences themselves if -c is set) paired with mates in <m1>
  <r>     Comma-separated list of files containing Crossbow-style reads.  Can be
          a mixture of paired and unpaired.  Specify "-" for stdin.
  <s>     Comma-separated list of files containing unpaired reads, or the
          sequences themselves, if -c is set.  Specify "-" for stdin.
  <hit>   File to write hits to (default: stdout)
Input:
  -q                 query input files are FASTQ .fq/.fastq (default)
  -f                 query input files are (multi-)FASTA .fa/.mfa
  -r                 query input files are raw one-sequence-per-line
  -c                 query sequences given on cmd line (as <mates>, <singles>)
  -C                 reads and index are in colorspace
  -Q/--quals <file>  QV file(s) corresponding to CSFASTA inputs; use with -f -C
  --Q1/--Q2 <file>   same as -Q, but for mate files 1 and 2 respectively
  -s/--skip <int>    skip the first <int> reads/pairs in the input
  -u/--qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)
  -5/--trim5 <int>   trim <int> bases from 5' (left) end of reads
  -3/--trim3 <int>   trim <int> bases from 3' (right) end of reads
  --phred33-quals    input quals are Phred+33 (default)
  --phred64-quals    input quals are Phred+64 (same as --solexa1.3-quals)
  --solexa-quals     input quals are from GA Pipeline ver. < 1.3
  --solexa1.3-quals  input quals are from GA Pipeline ver. >= 1.3
  --integer-quals    qualities are given as space-separated integers (not ASCII)
Alignment:
  -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities
    or
  -n/--seedmms <int> max mismatches in seed (can be 0-3, default: -n 2)
  -e/--maqerr <int>  max sum of mismatch quals across alignment for -n (def: 70)
  -l/--seedlen <int> seed length for -n (default: 28)
  --nomaqround       disable Maq-like quality rounding for -n (nearest 10 <= 30)
  -I/--minins <int>  minimum insert size for paired-end alignment (default: 0)
  -X/--maxins <int>  maximum insert size for paired-end alignment (default: 250)
  --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)
  --nofw/--norc      do not align to forward/reverse-complement reference strand
  --maxbts <int>     max # backtracks for -n 2/3 (default: 125, 800 for --best)
  --pairtries <int>  max # attempts to find mate for anchor hit (default: 100)
  -y/--tryhard       try hard to find valid alignments, at the expense of speed
  --chunkmbs <int>   max megabytes of RAM for best-first search frames (def: 64)
Reporting:
  -k <int>           report up to <int> good alignments per read (default: 1)
  -a/--all           report all alignments per read (much slower than low -k)
  -m <int>           suppress all alignments if > <int> exist (def: no limit)
  -M <int>           like -m, but reports 1 random hit (MAPQ=0); requires --best
  --best             hits guaranteed best stratum; ties broken by quality
  --strata           hits in sub-optimal strata aren't reported (requires --best)
Output:
  -t/--time          print wall-clock time taken by search phases
  -B/--offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)
  --quiet            print nothing but the alignments
  --refout           write alignments to files refXXXXX.map, 1 map per reference
  --refidx           refer to ref. seqs by 0-based index rather than name
  --al <fname>       write aligned reads/pairs to file(s) <fname>
  --un <fname>       write unaligned reads/pairs to file(s) <fname>
  --max <fname>      write reads/pairs over -m limit to file(s) <fname>
  --suppress <cols>  suppresses given columns (comma-delim'ed) in default output
  --fullref          write entire ref name (default: only up to 1st space)
Colorspace:
  --snpphred <int>   Phred penalty for SNP when decoding colorspace (def: 30)
     or
  --snpfrac <dec>    approx. fraction of SNP bases (e.g. 0.001); sets --snpphred
  --col-cseq         print aligned colorspace seqs as colors, not decoded bases
  --col-cqual        print original colorspace quals, not decoded quals
  --col-keepends     keep nucleotides at extreme ends of decoded alignment
SAM:
  -S/--sam           write hits in SAM format
  --mapq <int>       default mapping quality (MAPQ) to print for SAM alignments
  --sam-nohead       supppress header lines (starting with @) for SAM output
  --sam-nosq         supppress @SQ header lines for SAM output
  --sam-RG <text>    add <text> (usually "lab=value") to @RG line of SAM header
Performance:
  -o/--offrate <int> override offrate of index; must be >= index's offrate
  -p/--threads <int> number of alignment threads to launch (default: 1)
  --mm               use memory-mapped I/O for index; many 'bowtie's can share
  --shmem            use shared mem for index; many 'bowtie's can share
Other:
  --seed <int>       seed for random number generator
  --verbose          verbose output (for debugging)
  --version          print version information and quit
  -h/--help          print this usage message

