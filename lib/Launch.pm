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
use File::Temp qw/mktemp/;

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
            $logger->logdie("If placeholder <?> is used, exactly one expected file can be given");
        }
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
    
    if (0==system($cmd)){
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
