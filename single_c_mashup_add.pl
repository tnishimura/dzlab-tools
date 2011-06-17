#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use DBI;
use File::Basename;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if (! $opt_database || !$opt_gff);

if (!$opt_nickname){
    $opt_nickname = basename($opt_gff);
}

#######################################################################
# DBI

my $dbh = DBI->connect("dbi:SQLite:dbname=$opt_database","","", {RaiseError => 1, AutoCommit => 1});
my $check_sth = $dbh->prepare("select *  from scores where chr = ? and coord = ?");

$dbh->do("PRAGMA automatic_index = OFF");
$dbh->do("PRAGMA journal_mode = OFF");
$dbh->do("PRAGMA cache_size = 80000");

#$dbh->do("create table scores (chr, coord)");
#$dbh->do("create table columns (filename, nickname, num)");

#######################################################################
# Support

sub has_coord{
    my ($chr, $coord) = @_;
    $check_sth->execute($chr, $coord);
    return $check_sth->fetchrow_hashref;
}

sub colcount{
    my $dbh = shift;
    my ($count) = $dbh->selectrow_array("select count(*) from columns");
    return $count;
}

sub addcolumn{
    my ($dbh,$filename,$nickname) = @_;
    my $colnum = colcount($dbh)+1;
    $dbh->do("alter table scores add column _c_$colnum");
    $dbh->do("alter table scores add column _t_$colnum");
    $dbh->prepare("insert into columns (filename, nickname, num) values (?,?,?)")->execute($filename, $nickname, $colnum);
    return $colnum;
}

#######################################################################
# Main Loop


my $colnum = addcolumn($dbh,$opt_gff, $opt_nickname);
my $p = GFF::Parser->new(file => $opt_gff, normalize => 0);

my $insert_sth = $dbh->prepare("insert into scores(chr,coord,_c_$colnum,_t_$colnum) values (?,?,?,?)");
my $update_sth = $dbh->prepare("update scores set _c_$colnum = ?, _t_$colnum = ? where chr = ? and coord = ?");

my $counter = 0;

$dbh->{AutoCommit} = 0;
while (defined(my $gff = $p->next())){
    my ($chr, $coord, $c, $t) = ($gff->sequence,$gff->start,$gff->get_column('c'),$gff->get_column('t'));
    if ($colnum == 1){
        $insert_sth->execute($chr,$coord,$c,$t);
    }
    else{
        if (has_coord($chr, $coord)){
            $update_sth->execute($c,$t,$chr,$coord);
        }
        else{
            $insert_sth->execute($chr,$coord,$c,$t);
        }
    }
    if ($counter++ % 10000 == 0){
        $dbh->commit;
        if ($opt_verbose){
            say $counter;
        }
    }
    #warn $counter;
}
$dbh->commit;
$dbh->{AutoCommit} = 1;



$dbh->disconnect;

=head1 NAME
 
single_c_mashup_add.pl
 
=head1 SYNOPSIS

 single_c_mashup_add.pl

=head1 OPTIONS

=over

=item  -d <file> | --database <file>

=item  -g <single> | --gff <single>

=for Euclid
    single.type:        readable

=item  -n <nickname> | --nickname <nickname>

=item  -v | --verbose

=back

=cut

