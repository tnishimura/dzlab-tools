package Test::GFF;
use strict;
use Test::Builder::Module;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use GFF;

our @ISA    = qw(Test::Builder::Module);
our @EXPORT = qw(gff_has not_gff_has);
my $CLASS = __PACKAGE__;

sub _check{
    my $gff = shift;
    my %expected = %{shift()};

    my @errors;

    while (my ($key,$exp) = each %expected) {
        my $got = $gff->get_column($key);
        $got = uc $got if ($key eq 'sequence' and defined $got);
        $exp = uc $exp if ($key eq 'sequence' and defined $exp);
        if (defined $exp && ! defined $got){
            push @errors, "$key expected defined, got undefined\n";
        }
        elsif (! defined $exp && defined $got){
            push @errors, "$key expected undefined, got defined\n";
        }
        elsif (defined $exp && defined $got && $exp ne $got){
            push @errors, "$key expected $exp, got $got\n";
        }
    }
    return @errors;
}

sub gff_has{
    my $tb = $CLASS->builder;
    my $msg = pop;
    my @errors = _check(@_);
    #$tb->diag(@errors);
    return $tb->ok(0 == scalar(@errors), $msg) || $tb->diag(@errors);
}

sub not_gff_has{
    my $tb = $CLASS->builder;
    my $msg = pop;
    my @errors = _check(@_);
    return $tb->ok(0 != scalar(@errors), $msg) || $tb->diag(@errors);
}

1;

