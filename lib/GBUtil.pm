package GBUtil;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Bio::Graphics::Wiggle::Loader;
use File::Basename qw/basename/;
use Params::Validate qw/:all/;
use Hash::Util qw/lock_keys/;
use File::Which;
use Parallel::ForkManager;
use File::Path qw/make_path/;
use File::Spec::Functions qw/catfile rel2abs/;
use List::MoreUtils qw/notall/;
use YAML qw/LoadFile/;
use Readonly;

use FastaReader;
use GFF::Parser;
use GFF::Statistics qw/gff_detect_width/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw( load_mysql_config );

# my ($user, $pass, $database, $host) = load_mysql_config();
sub load_mysql_config{
    my $file = shift // "$ENV{HOME}/.bioperl";
    croak "no such file $file" unless -f $file;

    my $config = LoadFile($file);
    croak "can't read $file, malformatted?" unless ref $config eq 'HASH';
    # if (notall { exists $config->{$_} } qw/user pass database host/){
    #     croak "$file doesn't contain all of user, pass, database, host";
    # }
    return @{$config}{qw/user pass database host/};
}

1;
