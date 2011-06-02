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
use IPC::Cmd qw/run_forked/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(launch);

=head2 launch

 expected - This is a file (or arrayref of files) which we expect to be produce.
 force    - Run even if file exists.
 dryrun   - Don't actually run.

=cut

sub launch{
    my ($cmd, %opt) = @_;

    my $logger    = get_logger("PipeLine");

    if (!IPC::Cmd->can_use_run_forked()){
        $logger->logdie("can't run forked?");
    }

    my $force     = delete $opt{force} // 0;
    my $dryrun    = delete $opt{dryrun} // 0;

    $logger->info(join ", ", "running [$cmd]", ($force ? "forced" : ()), ($dryrun ? "dryrun" : ()));

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
            $tempfile = mktemp($expected[0] . '.tmpXXXXX');
            $cmd =~ s/\?\?/$tempfile/;
        }
        else {
            $logger->logdie("If placeholder ?? is used, exactly one expected file can be given");
        }
    }

    die "unknown parameters passed to doit" . Dumper \%opt if (%opt);

    if (!$force){
        if (! @expected){
            # none expected
        } elsif(@expected && grep {-f} @expected){
            $logger->info("Already done, skipping: [$cmd] ");
            return 1;
        }
    }
    if ($dryrun){
        $logger->info("Dryrun, exiting");
        return;
    }

    my $rv = run_forked($cmd, {
            discard_output => 1,
            stdout_handler => sub{ $logger->debug("stdout: " . shift); },
            stderr_handler => sub{ $logger->debug("stderr: " . shift); },
        });
    
    if ($rv->{exit_code} == 0){
        my $exp = join ", ", @expected;
        if (! @expected){
            $logger->info("Successfully launched and finished [$cmd]");
        } 
        elsif ($placeholder && -f $tempfile){
            rename $tempfile, $expected[0];
            $logger->info("Successfully launched and finished. Produced $expected[0] [$cmd]");
        } 
        elsif ($placeholder && ! -f $tempfile){
            $logger->logdie("command seems to have run but expected files $exp not produced [$cmd]");
        } 
        elsif (grep {-f} @expected){ 
            $logger->info("Successfully launched and finished. Produced $exp [$cmd]");
        } 
        else {
            $logger->logdie("command seems to have run but expected files $exp not produced [$cmd]");
        }
    } else {
        $logger->logdie("failed to run, dying: [$cmd]");
    }
}
1;
