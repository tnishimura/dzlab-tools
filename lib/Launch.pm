package Launch;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use Cwd;
use IPC::Open3;
use Log::Log4perl qw/get_logger/;
use Parallel::ForkManager;
use IPC::Cmd qw[can_run run ];

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(launch plaunch);

=head2 launch

 expected - This is a file (or arrayref of files) which we expect to be produce.
 force    - Run even if file exists.
 dryrun   - Don't actually run.

=cut

sub launch{
    my ($cmd, %opt) = @_;
    my $logger    = get_logger("PipeLine");

    my $force     = delete $opt{force} // 0;
    my $dryrun    = delete $opt{dryrun} // 0;

    $logger->info(join ", ", "running [$cmd]", ($force ? "forced" : ()), ($dryrun ? "dryrun" : ()));

    my $exe = (split ' ', $cmd)[0];
    if (!can_run($exe)){
        $logger->logdie("command not found: $exe");
    }

    ### commands can be arrayrefs or strings ###

    my @expected;
    if (exists $opt{expected}){
        if (ref $opt{expected} eq 'ARRAY'){ 
            @expected = @{$opt{expected}}
        } else {
            @expected = ($opt{expected});
        }
        delete $opt{expected}
    }

    die "unknown parameters passed to doit" . Dumper \%opt if (%opt);

    if (!$force){
        if (! @expected){
            # none expected
        } elsif(@expected && grep {-f} @expected){
            $logger->info("Already done, skipping: '$cmd' ");
            return 1;
        }
    }
    if ($dryrun){
        $logger->info("Dryrun, exiting");
        return;
    }

    ### in list context ###
    my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 1 );

    # a little bit dirty b/c the output will be written twice....
    $logger->info("output: " . join "", @$full_buf);

    if ($success){
        my $exp = join ", ", @expected;
        if (! @expected){
            $logger->info("Successfully launched and finished [$cmd]");
        } elsif (@expected == grep {-f} @expected){ 
            $logger->info("Successfully launched and finished. Produced $exp [$cmd].");
        } else {
            # mark all expected files that were produced as dirty
            for my $e (@expected) { if (-f $e){ rename $e, "DIRTY_$e"; } }
            $logger->logdie("command seems to have run but expected files $exp not produced [$cmd]");
        }
    } else {
        if (@expected){
            # mark all expected files that were produced as dirty
            for my $e (@expected) { if (-f $e){ rename $e, "DIRTY_$e"; } }
        }
        $logger->logdie("failed to run, dying: [$cmd]");
    }
}
sub plaunch{
    my ($numprocs, @jobs) = @_;
    # Max 30 processes for parallel download
    my $pm = new Parallel::ForkManager($numprocs);

    foreach my $j (@jobs) {
        $pm->start and next; # do the fork
        launch(@$j);
        $pm->finish; # do the exit in the child process
    }
    $pm->wait_all_children;
}
1;
