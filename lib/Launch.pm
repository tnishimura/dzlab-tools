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
use IPC::Cmd qw/run/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(launch);

=head2 launch

 expected - This is a file (or arrayref of files) which we expect to be produce.
 force    - Run even if file exists.
 dryrun   - Don't actually run.
 also     - Also print output of command to this file.

=cut

sub launch{
    my ($cmd, %opt) = @_;

    my $logger = get_logger();

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
            $logger->logdie("If placeholder ?? is used, exactly one expected file can be given");
        }
    }

    $logger->info(join ", ", "running [$cmd]", ($force ? "forced" : ()), ($dryrun ? "dryrun" : ()));

    die "unknown parameters passed to launch" . Dumper \%opt if (%opt);

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

    $logger->info("Running $cmd");
    my ($success, $errmsg, $fullbuf) = run(command => $cmd, verbose => 1);

    if (! $success){
        $logger->logwarn("FAILED: $cmd");
        $logger->logwarn("$errmsg");
        $logger->logdie("dying...");
    }

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
