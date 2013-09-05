package PFMLauncher;
use v5.12.0;
use MooseX::Singleton;
use Data::Dumper;
use Carp;
use autodie;    
use Parallel::ForkManager;
use File::Temp qw/mktemp/;
use POSIX qw/strftime/;
use File::Copy qw/move/;
use Sys::Info;

has processes => (
    is => 'ro',
    required => 0, 
    default => sub { _cpucount() },
);

has pfm => (
    is => 'ro', 
    lazy_build => 1,
);

has failed_cmd_pids => (
    traits  => ['Array'],
    is      => 'bare',
    isa     => 'ArrayRef[Int]',
    default => sub { [] },
    handles => {
        failed_cmd_pids    => 'elements',
        add_failed_cmd_pid => 'push',
    },
);

has pid_to_cmd => (
    traits    => ['Hash'],
    is        => 'ro',
    isa       => 'HashRef[Str]',
    default   => sub { {} },
    handles   => {
        set_pid_to_cmd => 'set',
        get_pid_to_cmd => 'get',
    },
);


sub _build_pfm{
    my $self = shift;
    my $pfm = Parallel::ForkManager->new($self->processes);

    $pfm->run_on_finish(sub{
            my ($pid, $exit_code, $ident) = @_;
            if ($exit_code != 0){
                $self->add_failed_cmd_pid($pid);
            }
        });
    return $pfm;
}

sub wait_all_children{
    my $self = shift;
    $self->pfm->wait_all_children;
    my @failed_cmd = map { $self->get_pid_to_cmd($_) } $self->failed_cmd_pids;
    if (@failed_cmd){
        _logdie("following commands failed:", @failed_cmd);
    }
}

sub launch_and_wait{
    my ($self, $cmd, %opt) = @_;
    $opt{wait} = 1;
    $self->launch($cmd, %opt);
}

sub launch{
    my ($self, $cmd, %opt) = @_;

    my $wait        = delete $opt{wait};
    my $stdin_spec  = delete $opt{stdin};
    my $stdout_spec = delete $opt{stdout};
    my $stderr_spec = delete $opt{stderr};
    my @expected;
    if (exists $opt{expected}){
        if (ref $opt{expected} eq 'ARRAY'){
            @expected = @{$opt{expected}};
        }
        elsif (! ref $opt{expected}){
            @expected = ($opt{expected});
        }
    }
    delete $opt{expected};

    use List::MoreUtils qw/all/;
    
    if (@expected && all { -f } @expected){
        _info("All expected files listed below already existed, not rerunning", @expected);
        $self->wait_all_children if $wait;
        return;
    }

    my $placeholder = $cmd =~ /%/;
    my $tempfile;
    if ($placeholder){
        if (@expected == 1){
            $tempfile = mktemp($expected[0] . '.tmp.XXXXX');
            $cmd =~ s/%/"$tempfile"/;
        }
        else {
            _logdie("If placeholder % is used, exactly one expected file can be given");
        }
    }

    die "unknown parameters passed to launch" . Dumper \%opt if (%opt);

    my $pfm = $self->pfm;
    my $pid = $pfm->start; 
    
    if ($pid != 0){ # in parent
        $self->set_pid_to_cmd($pid, $cmd);
        $self->wait_all_children if $wait;
        return;
    }
    else{ # in child
        _setup_input_fh (\*STDIN , $stdin_spec  );
        _setup_output_fh(\*STDOUT, $stdout_spec );
        _setup_output_fh(\*STDERR, $stderr_spec );

        _info("launching [$cmd]");

        my $rc = system($cmd);

        if ($rc == 0){
            my $exp = join ", ", @expected;
            if (! @expected){
                _info("Successfully launched and finished [$cmd]");
                $pfm->finish(0);
            } 
            elsif ($placeholder && -f $tempfile){
                if (move $tempfile, $expected[0]){
                    _info("Successfully launched and finished. Produced $expected[0] [$cmd]");
                    $pfm->finish(0);
                }
                else{
                    _info("Successfully launched and finishedbut couldn't rename $tempfile to $expected[0]? [$cmd]");
                    $pfm->finish(1);
                }

            } 
            elsif ($placeholder && ! -f $tempfile){
                _logdie("command seems to have run but expected files $exp not produced [$cmd]");
                $pfm->finish(1);
            } 
            elsif (@expected == grep {-f} @expected){ 
                _info("Successfully launched and finished. Produced $exp [$cmd]");
                $pfm->finish(0);
            } 
            else {
                _logdie("command seems to have run but expected files $exp not produced [$cmd]");
                $pfm->finish(1);
            }
        } else {
            _logdie("failed to run, dying: [$cmd]");
            $pfm->finish($rc);
        }
    }

    
    return 1;
}


sub _setup_input_fh{
    my $target_fh = shift;
    my $fh_or_filename = shift;
    if (! defined $fh_or_filename){
        return;
    }
    else{
        close $target_fh;
        if (ref $fh_or_filename eq 'GLOB'){
            $target_fh = $fh_or_filename;
        }
        else {
            $fh_or_filename =~ s/^<//;
            open $target_fh, '<', $fh_or_filename;
        }
    }
}

sub _setup_output_fh{
    my $target_fh = shift;
    my $fh_or_filename = shift;
    if (! defined $fh_or_filename){
        return;
    }
    else{
        close $target_fh;
        if (ref $fh_or_filename eq 'GLOB'){
            $target_fh = $fh_or_filename;
        }
        elsif ($fh_or_filename =~ s/^>>//){
            open $target_fh, '>>', $fh_or_filename;
        }
        else{
            $fh_or_filename =~ s/^>//; 
            open $target_fh, '>', $fh_or_filename;
        }
    }
}

sub _cpucount{
    my $info = Sys::Info->new;
    my $cpu = $info->device('CPU');
    return $cpu->count;
}

# add real logging
sub _logdie{
    my @msgs = @_;
    my $time = strftime("%H:%M:%S", localtime(time()));
    for my $msg (@msgs) { STDERR->print("[$time] $msg\n"); }
    exit(1);
}
sub _info {
    my @msgs = @_;
    my $time = strftime("%H:%M:%S", localtime(time()));
    for my $msg (@msgs) { STDERR->print("[$time] $msg\n"); }
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;


__END__

=head1 PFMLauncher.pl 

 use aliased PFMLauncher => 'PFM';

 PFM->launch_and_wait("ls -l"); # adds wait => 1 to opt

 PFM->launch("sleep 1; ls -l", stderr => 'qwer', stdout => 'asdf');
 PFM->launch("sleep 1; echo foo > %", expected => ['asdf2'], wait => 1);
 PFM->launch("sleep 1; asdfasdf");
 PFM->launch("sleep 1; touch a; touch b", expected => ['a', 'b']);
 PFM->launch("sleep 1; touch a; touch b", expected => ['a', 'b','c']);
 PFM->wait_all_children();

=cut

