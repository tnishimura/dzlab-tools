package Memo;
use version; 
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(pmemo);

use CHI;
use Digest::file qw//;
use YAML qw//;
use English;

our $cache = CHI->new( driver => 'File',
    root_dir => '/tmp/DZ'
);

sub pmemo{
    my ($filename, $propertyname, $sub) = @_;

    my $md5 = Digest::file::digest_file_hex($filename, "MD5");
    my $key = $md5 . $SUBSCRIPT_SEPARATOR . $propertyname;

    if (defined (my $value = $cache->get($key))){
        return YAML::Load $value;
    }
    else{
        local $_ = $filename;
        my $value = $sub->();
        $cache->set($key, YAML::Dump $value);
        return $value;
    }
}
# say Dumper memop "/etc/passwd", "size5", sub{ [123, 456, { a => 1, x => 312 } ]};

1;
