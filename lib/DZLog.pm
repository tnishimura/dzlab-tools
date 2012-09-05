package DZLog;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Log::Dispatch;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(initialize_dzlog INFO WARN DEBUG ERROR);

my $logger;
my @default_args = ( outputs => [ [ 'Screen', min_level => 'warning' ], ] );

sub initialize_dzlog{
    if (! $logger){
        $logger = Log::Dispatch->new(0 == @_ ? @default_args : @_);
    }
    else {
        carp "double initialization of " . __PACKAGE__ ;
    }
}


sub INFO  { initialize_dzlog() unless $logger; $logger->info(@_) };
sub WARN  { initialize_dzlog() unless $logger; $logger->warn(@_) };
sub DEBUG { initialize_dzlog() unless $logger; $logger->debug(@_) };
sub ERROR { initialize_dzlog() unless $logger; $logger->error(@_); exit 1};

1;
