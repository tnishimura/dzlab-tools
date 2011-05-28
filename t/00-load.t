#!perl -T

use Test::More tests => 1;

BEGIN {
    use_ok( 'Tree::Range' ) || print "Bail out!\n";
}

diag( "Testing Tree::Range $Tree::Range::VERSION, Perl $], $^X" );
