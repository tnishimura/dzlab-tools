package Launch;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use Cwd;
use Log::Log4perl qw/get_logger/;
use Parallel::ForkManager;
use File::Temp qw/mktemp/;
use autodie;
use IPC::Cmd qw/run/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(launch);

our $VERBOSE = exists $ENV{QUIET} ? ! $ENV{QUIET} : 1;

=head2 launch

 expected - This is a file (or arrayref of files) which we expect to be produce.
 force    - Run even if file exists.
 dryrun   - Don't actually run.
 also     - Also print output of command to this file.

=cut

{
    sub _logdie{
        if (Log::Log4perl->initialized()){ 
            my $logger = get_logger("PipeLine");
            $logger->logdie($_[0]); 
        }
        else{ 
            die $_[0]; 
        }
    }
    sub _info{
        if (Log::Log4perl->initialized()){ 
            my $logger = get_logger("PipeLine");
            $logger->info($_[0]); 
        }
        elsif ($VERBOSE){
            say STDERR $_[0]; 
        }
    }
}

sub launch{
    my ($cmd, %opt) = @_;

    my $force     = delete $opt{force} // 0;
    my $dryrun    = delete $opt{dryrun} // 0;
    my $also      = delete $opt{also} // 0;

    my @expected;
    if (exists $opt{expected}){
        if (ref $opt{expected} eq 'ARRAY'){ 
            @expected = @{$opt{expected}}
        } else {
            @expected = ($opt{expected});
        }
        delete $opt{expected}
    }

    my $placeholder = $cmd =~ /\?\?/;
    my $tempfile;
    if ($placeholder){
        if (@expected == 1){
            $tempfile = mktemp($expected[0] . '.tmp.XXXXX');
            $cmd =~ s/\?\?/$tempfile/;
        }
        else {
            _logdie("If placeholder ?? is used, exactly one expected file can be given");
        }
    }


    die "unknown parameters passed to launch" . Dumper \%opt if (%opt);

    if (!$force){
        if (! @expected){
            # none expected
        } elsif(@expected && grep {-f} @expected){
            _info("Already done, skipping: [$cmd] ");
            return 1;
        }
    }
    if ($dryrun){
        _info("Dryrun, exiting");
        return;
    }

    # no need to say _info this b/c verbose => 1 does it for us.
    my ($success, $errmsg, $fullbuf) = run(command => $cmd, verbose => $VERBOSE && IPC::Cmd->can_use_ipc_run());

    if (! $success){
        _logdie("FAILED: $cmd\n$errmsg\ndying...");
    }

    _info(join "", @$fullbuf);

    if (@$fullbuf && $also){
        open my $also_fh, '>>', $also;

        print $also_fh "===== Output of [$cmd] was:\n";
        print $also_fh join "", @$fullbuf;
        print $also_fh "=====\n";

        close $also_fh if $also;
    }

    if ($success == 1){
        my $exp = join ", ", @expected;
        if (! @expected){
            _info("Successfully launched and finished [$cmd]");
        } 
        elsif ($placeholder && -f $tempfile){
            rename $tempfile, $expected[0];
            _info("Successfully launched and finished. Produced $expected[0] [$cmd]");
        } 
        elsif ($placeholder && ! -f $tempfile){
            _logdie("command seems to have run but expected files $exp not produced [$cmd]");
        } 
        elsif (grep {-f} @expected){ 
            _info("Successfully launched and finished. Produced $exp [$cmd]");
        } 
        else {
            _logdie("command seems to have run but expected files $exp not produced [$cmd]");
        }
    } else {
        _logdie("failed to run, dying: [$cmd]");
    }
    return 1;
}
1;
