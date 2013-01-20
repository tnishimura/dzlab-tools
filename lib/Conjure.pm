package Conjure;
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Carp;
use autodie;    
use Data::Dumper;
use POE; 
use POE::Wheel::Run; 
use IO::Handle;
use Params::Validate qw//;
use Mouse;
use Mouse::Exporter;
 
Mouse::Exporter->setup_import_methods( as_is => [ 'conjure' ],);

our $DEBUG = 0;
sub debug{ say $_[0] if $DEBUG; }

has program        => (is => 'rw', required => 1, isa => 'Str|ArrayRef[Str]|CodeRef');
has on_stderr      => (is => 'ro', required => 0, isa => 'CodeRef');
has on_stdout      => (is => 'ro', required => 0, isa => 'CodeRef');
has stdin          => (is => 'ro', required => 0, isa => 'GlobRef');
has pid            => (is => 'rw', init_arg => undef);
has exit_status    => (is => 'rw', init_arg => undef);

# Private to conjure() function
has tee_stdout_fh => (is => 'ro', required => 0, isa => 'GlobRef');
has tee_stderr_fh => (is => 'ro', required => 0, isa => 'GlobRef');

sub BUILD{
    my ($self) = @_;

    POE::Session->create(
        object_states => [
            $self => [qw{ _start _stderr_handler _stdout_handler _stop _close_event _sig_child}],
        ],
    );

    return $self;
}

sub _start{ 
    my ($kernel, $heap, $session, $self) = @_[KERNEL, HEAP, SESSION, OBJECT];

    my $wheel = POE::Wheel::Run->new(
        Program       => $self->program,
        StdoutEvent   => '_stdout_handler',
        StderrEvent   => '_stderr_handler',
        CloseEvent    => '_close_event',
        ($self->stdin ? (RedirectStdin => $self->stdin) : ()),
    );

    my $wid = $wheel->ID;
    debug "making a run wheel $wid";

    $kernel->sig_child($wheel->PID, '_sig_child');
    $heap->{wheel} = $wheel;
    $self->pid($wheel->PID);
}

sub _stop{
    debug "_stop";
    my ($heap, $self) = ($_[HEAP], $_[OBJECT]);
} 

# separate from _close_event b/c _sig_child gets passed 
sub _sig_child{
    my ($heap, $self, $pid, $exit_status) = ($_[HEAP], $_[OBJECT], $_[ARG1], $_[ARG2]);
    debug "_sig_child pid $_[ARG1] exited with status $_[ARG2].\n";
    $self->pid($pid);
    $self->exit_status($exit_status);

    delete $heap->{wheel};
}

sub _close_event {
    my ($heap, $self) = ($_[HEAP], $_[OBJECT]);
    debug "_close_event";

    delete $heap->{wheel};
}

our $PID;
sub _stdout_handler {
    my ($line, $self) = @_[ARG0, OBJECT];
    debug "_stdout_handler got: $line";

    if ($self->on_stdout){
        local $_ = $line;
        local $PID = $self->pid;
        $self->on_stdout->();
    }
    if ($self->tee_stdout_fh){
        $self->tee_stdout_fh->say($line);
    }
}

sub _stderr_handler {
    my ($line, $self) = @_[ARG0, OBJECT];
    debug "_stderr_handler got: $line";

    if ($self->on_stderr){
        local $_ = $line;
        local $PID = $self->pid;
        $self->on_stderr->();
    }
    if ($self->tee_stderr_fh){
        $self->tee_stderr_fh->say($line);
    }
}

=head2 conjure()

 my ($ok, $pid, $exit_status, $msg) = 
 conjure(
     program => ['/bin/ls', '-l'], 
     on_std_out => sub { ... },
     on_std_err => sub { ... },
     stdin => FILEHANDLE
 );

 my $ok = conjure(
     program => ['/bin/ls', '-l'], 
     on_std_out => sub { ... },
     on_std_err => sub { ... },
     stdin => FILEHANDLE
 );

=cut
sub conjure{
    my %opt = Params::Validate::validate(@_, {
            tee_stdout => 0,

            # four arguments passed directly to Conjure->new();
            program   => 1,
            on_stderr => 0,
            on_stdout => 0,
            stdin     => 0,
        });

    my $program_string = _program_string($opt{program});

    #######################################################################
    # prepare tee file handles

    my $tee_stdout = delete $opt{tee_stdout};
    
    # check if file exists
    if (defined $tee_stdout && -f $tee_stdout && -s $tee_stdout > 0){
        if (wantarray){
            return (1, -1, 0, "[$program_string]: $tee_stdout already exists and is non-zero in size, not conjuring");
        }
        else{
            return 1;
        }
    }

    # create tmp file name and file handle
    my ($tee_stdout_tmp, $tee_stdout_fh) = do{
        if (defined $tee_stdout){
            my $tmp = "$tee_stdout.tmp" . _randstring(6);
            open my $fh, '>', $tmp;
            ($tmp, $fh);
        }
        else{
            (undef, undef);
        }
    };

    # say Dumper {
    #     tee_stdout     => $tee_stdout,
    #     tee_stdout_tmp => $tee_stdout_tmp,
    #     tee_stdout_fh  => $tee_stdout_fh,
    # };

    #######################################################################
    # run
    
    # according to POE::Wheel::Run doc, this is necessary to allow nesting
    POE::Kernel->stop(); 

    my $conjure = Conjure->new(
        %opt,
        (defined $tee_stdout_fh ? (tee_stdout_fh => $tee_stdout_fh) : ())
    );
    POE::Kernel->run();

    #######################################################################
    # exit
    
    my $exit = $conjure->exit_status();
    my $pid = $conjure->pid();
    my $nick = "[$program_string ($pid)]";

    # always close tee fh
    close $tee_stdout_fh if (defined $tee_stdout_fh);

    # child didn't exit correctly
    if ($exit != 0){
        return wantarray ?  (0, $pid, $exit, "$nick exited with exit status $exit") : return 0;
    }
    # exitted correctly
    else{
        # expecting tee'd file to exist:
        if (defined $tee_stdout_fh){
            # and it exists properly:
            if (-f $tee_stdout_tmp && -s $tee_stdout_tmp > 0){
                rename $tee_stdout_tmp, $tee_stdout;
                return wantarray ? (1, $pid, $exit, "$nick exited successfully, $tee_stdout created successfully"): 1;
            }
            # but file wasn't created properly?
            else{
                return wantarray ? (0, $pid, $exit, "$nick exited successfully, but $tee_stdout was not created?"): 0;
            }
        }
        # no tee'ing, successful exit code
        else{
            return wantarray ? (1, $pid, $exit, "$nick successfully"): 1;
        }
    }
}

# random a-z chars. halfassed but better than bringing in new module
sub _randstring{
    my $num = shift;
    join "", map {
        chr(int(rand(26)) + 97)
    } (1 .. $num);
}

sub _program_string{
    my $prog = shift;
    if (ref $prog eq 'ARRAY'){
        return join " ", @$prog;
    }
    elsif (ref $prog eq 'CODE'){
        return "sub { ... }";
    }
    else{
        return $prog;
    }
}

1;
