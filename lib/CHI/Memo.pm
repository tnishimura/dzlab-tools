package CHI::Memo;
use version; 
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(dump clear_cache);
our @EXPORT = qw(chimemo);

use CHI;
use Digest::file qw//;
use YAML qw//;
use English;

our $root_dir = '/tmp/DZ';
our $cache; 

sub clear_cache{
    if (! defined $cache){ _init(); }
    $cache->clear();
}
sub dump{
    if (! defined $cache){ _init(); }

    for my $k ($cache->get_keys()) {
        my ($md5, $property) = split /$SUBSCRIPT_SEPARATOR/, $k;
        say Dumper {"$k" => YAML::Load $cache->get($k)};
    }
}

sub _init{
    $cache = CHI->new( driver => 'File', root_dir => $root_dir,);
}

=head2 chimemo FILENAME, PROPERTY, SUB

chimemo('TAIR_reference.fas', 'lengths', sub { $fr->get_lengths() }

=cut
sub chimemo{
    my ($filename, $propertyname, $sub) = @_;
    if (! defined $cache){ _init(); }

    my $md5 = Digest::file::digest_file_hex($filename, "MD5");
    my $key = $md5 . $SUBSCRIPT_SEPARATOR . $propertyname;

    if (defined (my $value = $cache->get($key))){
        return @{YAML::Load $value};
    }
    else{
        local $_ = $filename;
        my @value = $sub->();
        $cache->set($key, YAML::Dump \@value);
        return @value;
    }
}
# say Dumper memop "/etc/passwd", "size5", sub{ [123, 456, { a => 1, x => 312 } ]};

1;
