#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use DBI;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use DZUtil qw/timestamp/;
use File::Temp qw/mktemp/;
use FastaReader;
use Counter;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) if $opt_help;

if (! $opt_no_skip && $opt_reference){
    warn "why are you giving me a reference file if you're not going to do --no-skip?";
}

my $counter = Counter->new();

#######################################################################
# Get lengths

my %lengths;
if ($opt_no_skip && $opt_reference){
    # normalize => 1 means make sure to uc all keys.... need to fix this
    my $fr = FastaReader->new(file => $opt_reference, normalize => 1); 
    %lengths = %{$fr->length};
}
if ($opt_verbose){
    say Dumper \%lengths;
}

#######################################################################
# Output

my $tmpfile;
if ($opt_output eq '-' && @opt_files == 1){
    my $defout = $opt_files[0] . ".merged";
    $tmpfile = mktemp("$defout.tmp.XXXXX");
    open STDOUT, '>', $defout or die "can't open $opt_output for writing";
} 
elsif ($opt_output eq '-' && @opt_files > 1){
    die "Sorry, when there are more than 1 input file, you need to specify a name with -o";
}
else {
    open my $fh, '>', $opt_output or die "can't open $opt_output for writing";
    $tmpfile = mktemp("$opt_output.tmp.XXXXX");
    select $fh;
}

#######################################################################
# Create Table

my $dbh = DBI->connect("dbi:SQLite:dbname=$tmpfile","","", {RaiseError => 1, AutoCommit => 0});
$dbh->do("PRAGMA automatic_index = OFF");
$dbh->do("PRAGMA journal_mode = OFF");
$dbh->do("PRAGMA cache_size = 80000");

$dbh->do(q{
    create table gff (sequence, position integer, c integer, t integer, n integer)
    });

$dbh->do("create index idx1 on gff (sequence,position)");
my $insert_sth = $dbh->prepare("insert into gff (sequence, position, c, t, n) values (?,?,?,?,1)");
my $update_sth = $dbh->prepare("update gff set c=?, t=?, n=? where sequence=? and position=?");
my $checker_sth = $dbh->prepare("select count(*) as present, c, t, n from gff where sequence=? and position=?");
#my %inserted = ();
sub record{
    my ($sequence, $position, $new_c, $new_t) = @_;
    my $key = $sequence . "~" . $position;

    #if (exists $inserted{$key}){

    #}
    
    $checker_sth->execute($sequence, $position); 
    my ($exists, $c, $t, $n) = $checker_sth->fetchrow_array;
    if ($exists){
        $update_sth->execute($c+$new_c, $t+$new_t, $n+1, $sequence, $position);
    }
    else{
        $insert_sth->execute($sequence, $position, $new_c, $new_t);
    }
}


say STDERR "tmp DB $tmpfile created";

#######################################################################
# Insert

if ($opt_verbose){
    say STDERR "started inserting at: " . timestamp();
}

my $feature;

#my $counter = 0;
#my $commit_size = 20000;
for my $file (@opt_files) {
    my $p = GFF::Parser->new(file => $file, normalize => 0);
    while (defined(my $gff = $p->next())){
        if (my $count = $counter->increment()){
            $dbh->commit;

        }
        if (defined $feature && defined $gff->feature && $feature ne $gff->feature){
            die "merging more than one feature at once? not supported yet";
        }

        $feature //= $gff->feature;

        my ($start, $end) = ($gff->start, $gff->end);
        if ($start != $end){
            die "$0 is only for single-c file";
        }
        # round up to the nearest window (101-150 to 150, 151-200 to 200, etc for w50)
        my $windowed_position = ($start - 1) + ($opt_window_size - ($start - 1) % $opt_window_size);
        record($gff->sequence, $windowed_position, ($gff->get_column('c')//0), ($gff->get_column('t')//0));
        #$insert_sth->execute($gff->sequence, $windowed_position, ($gff->get_column('c')//0), ($gff->get_column('t')//0));
    }
}

$dbh->commit;
$dbh->{AutoCommit} = 1;

if ($opt_verbose){
    say STDERR "done inserting at: " . timestamp();
}

#######################################################################
# Select

my $select = $dbh->prepare(<<SELECT );
    select sequence, position, c, t, n from gff 
SELECT

#my $select = $dbh->prepare(<<SELECT );
#    select sequence, position, sum(c) as c, sum(t) as t, count(*) as n
#    from gff 
#    group by sequence, position
#    order by sequence, position
#SELECT

$select->execute();

my %last_ends; 

while (defined(my $row = $select->fetchrow_hashref())){
    my ($sequence, $position, $c, $t, $n) = @{$row}{qw/sequence position c t n/};
    my ($current_start, $current_end) = ($position - $opt_window_size + 1, $position);
    my $score = $opt_report_count ? $n : sprintf "%.4f", ($c+$t == 0)? 0 : $c/($c+$t);

    if (!exists $last_ends{$sequence}){
        $last_ends{$sequence} = $opt_window_size;
    }

    # catch up
    if ($opt_no_skip){
        while ($last_ends{$sequence} < $position ){
            say join "\t", $sequence, '.', $feature, $last_ends{$sequence}-$opt_window_size+1, 
            $last_ends{$sequence}, ($opt_report_count ? 0 : '.'), ('.') x 3;
            $last_ends{$sequence}+=$opt_window_size;
        }
    }

    # if it's the last window, trim appropriately
    if (exists $lengths{uc $sequence} && $current_end >= $lengths{uc $sequence}){
        say join "\t", $sequence, '.', $feature, $current_start,$lengths{uc $sequence}, $score, '.', '.', "c=$c;t=$t;n=$n";
        $last_ends{$sequence} = $lengths{uc $sequence};
    }
    else{
        say join "\t", $sequence, '.', $feature, $current_start,$current_end, $score, '.', '.', "c=$c;t=$t;n=$n";
        $last_ends{$sequence}   = $current_end   + $opt_window_size;
    }
}

#######################################################################
# for all sequence that were seen and go up to the end

# bug: window 1 with no-skip and reference leads to final position not being printed? 

#say STDERR Dumper \%last_ends;
#say STDERR Dumper \%lengths;

if ($opt_no_skip && $opt_reference){

    # but only for unfinished seqs:
    for my $seq_seen (grep { $last_ends{$_} < $lengths{uc $_} } keys %last_ends) {

        while ($last_ends{$seq_seen} < $lengths{uc $seq_seen} ){
            say join "\t", $seq_seen, '.', $feature, $last_ends{$seq_seen}-$opt_window_size+1, $last_ends{$seq_seen}, ('.') x 4;
            $last_ends{$seq_seen}+=$opt_window_size;
        }

        # if length was not multiple of opt_window_size, there'd be one left over--
        my $l = $lengths{uc $seq_seen};
        if ($l % $opt_window_size != 0){
            say join "\t", $seq_seen, '.', $feature, 
            ($l - $l % $opt_window_size + 1),
            $l,
            ('.') x 4;
        }
    }
}

if ($opt_verbose){
    say STDERR "done at: " . timestamp();
}

$dbh->disconnect;
unlink $tmpfile;

=head1 NAME
 
window_by_fixed.pl
 
=head1 SYNOPSIS

 window_by_fixed.pl

=head1 OPTIONS

=over

=item  -w <size> | --window-size <size>

=for Euclid
    size.default:     1
    size.type:        int, size >= 1 
    size.type.error:  <size> must be integer greater than 1

=item  -r <fasta> | --reference <fasta>

Reference fasta file.  If --no-skip/-k is used, for empty windows between the last window and the end of the
sequences will also be outputted.

=for Euclid
    fasta.type:        readable

=item -o <file> | --output <file>

=for Euclid 
    file.default: '-'

=item <files>...

=item  -k | --no-skip

Don't omit windows without any scores.  Currently only works for files with single sequences.

=item -h | --help

=item -v | --verbose

=item  -n | --report-count 

Report 'n' in the scores column instead of c/(c+t).

=back

=cut

