use 5.006;
use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Tree::Range',
    license             => 'perl',
    dist_author         => q{<Author> <<email@localhost>>},
    dist_version_from   => 'lib/Tree/Range.pm',
    build_requires => {
        'Test::More' => 0,
        'perl' => 5.010,
        'Moose' => 0,
    },
    add_to_cleanup      => [ 'Tree-Range-*' ],
    create_makefile_pl => 'small',
);

$builder->create_build_script();
