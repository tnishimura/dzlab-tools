package Log::Flogger;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Carp qw/confess carp croak/;
use Fcntl qw/:flock SEEK_END/;
use File::Basename qw/basename/;
use Hash::Util qw/lock_keys/;
use POSIX qw/strftime/;
use YAML qw/Load Dump/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(debug info);

my %name2numeric = (DEBUG => 3,   INFO => 2,   WARN => 1,   FATAL => 0);
#my %numeric2name = (3 => 'DEBUG', 2 => 'INFO', 1 => 'WARN', 0 => 'FATAL');
my $scriptname = basename($0);

#######################################################################
# load configs or initialize

our %config;
if (defined $ENV{FLOGGER}){
    %config = %{Load $ENV{FLOGGER}};
}
else{
    %config = (files => [], level => 'INFO', starttime => time());
    $ENV{FLOGGER} = Dump \%config;
}
lock_keys %config;

#######################################################################
# set_file set_level

sub set_level{
    my $level = shift;
    confess "BUG" if (!defined $level);
    confess "invalid level $level" if (!exists $name2numeric{$level});

    $config{level} = $level;
    $ENV{FLOGGER} = Dump \%config;
}

sub set_file{
    my ($file) = @_;
    confess "BUG: sqlogger_init() file arg" if ! defined $file;

    push @{$config{files}}, $file;
    $ENV{FLOGGER} = Dump \%config;

    # every new file should have a starting timestamp
    _log('INFO', 1, strftime("starting log for $scriptname at %Y/%m/%d %H:%M:%S", localtime(time())));
}

#######################################################################
# _log

sub debug{ _log('DEBUG', 1, @_); }
sub info { _log('INFO', 1, @_); }
$SIG{__WARN__} = sub{ _sig_log('WARN', 2, @_); };
$SIG{__DIE__}  = sub{ _sig_log('FATAL', 2, @_); };

sub _sig_log{
    my ($level, $depth, $msg) = @_;
    chomp $msg; 
    $msg =~ s{at \Q$0\E line \d+\.$}{};
    _log($level, $depth, $msg);
}

sub _log{
    #confess "sqlog: sqlogger not initialized" if (! defined $_insert_sth);
    my ($level, $depth, $message) = @_;
    if (! defined $level || ! defined $message){
        confess "sqlog: BUG, usage error";
    }

    # always print if breaks threshold
    if ($name2numeric{$level} <= $name2numeric{$config{level}}){
        say STDERR "$level: $message";
    }

    # get the file and the line of the log_* function caller
    my ($file, $line) = (caller($depth))[1,2];
    $file = basename $file;

    # get the name of the function that called the log_*
    my $function;
    if (my @c = caller($depth + 1)){
        $function = $c[3];
    }
    else{
        $function = "__TOP__";
    }

    #strftime("%H:%M:%S", localtime(time()))
    my $logmsg = sprintf(
        "[%s - %s (%s) %s(%s) %s]\n%s\n", 
        $scriptname,
        _minutes(time() - $config{starttime}), 
        $level, 
        $file, 
        $line, 
        $function, 
        join("\n", map { "  $_" } split /\n/, $message),
    ); 

    for my $file (@{$config{files}}) {
        open my $fh, '>>', $file;

        # lock... don't die/warn here so we don't loop
        if (! flock($fh, LOCK_EX) || ! seek($fh, 0, SEEK_END)){
            say STDERR "Can't lock $file? Potential bug. quitting";
            exit 1;
        }

        print $fh $logmsg;

        # unlock
        if (! flock($fh, LOCK_UN)){
            say STDERR "Can't unlock $file? Potential bug. quitting";
            exit 1;
        }

        close $fh;
    }
}


sub _minutes{
    my $seconds = shift;
    my $m = int($seconds / 60);
    my $s = $seconds % 60;

    return sprintf("%d:%02d", $m, $s);
}

1; # End of Log::SQLogger
