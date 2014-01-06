#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use File::Find;
use File::Basename;
use File::Spec::Functions;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) if !$opt_directory;

if ($opt_output ne '-'){
    open my $fh, '>', $opt_output;
    select $fh;
}

my @freqs;
find( sub {
        # $File::Find::name - filename relative to pwd
        # $File::Find::dir  - dirname relative to pwd 
        # $_                - filename relative to $File::Find::dir
        
        if (basename($File::Find::dir) =~ /^single-c/ && m/freq$/){
            push @freqs, $File::Find::name;
        }
    }, $opt_directory);

#say STDERR Dumper \@freqs;

my $first = 1;

for my $f (sort @freqs) {
    #say STDERR $f;
    open my $fh, '<', $f;
    my ($header, $body) = <$fh>;
    chomp $header;
    chomp $body;
    close $fh;

    if ($first){
        $first = 0;
        say "Filename\t$header";
    } 
    say "$f\t$body";
}



=head1 NAME

collect-freqs.pl - ...

=head1 SYNOPSIS

Usage examples:

 collect-freqs.pl [options]...

=head1 OPTIONS

=over

=item  -o <freq> | --output <freq>

=for Euclid
    freq.default:     '-'

=item <directory>

=back

=cut

