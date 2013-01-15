package Sam::Parser;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Hash::Util qw/lock_keys unlock_keys/;
use List::MoreUtils qw/all/;
use Any::Moose;

extends 'Parser';

# from headers
has program_name    => (is => 'rw', init_arg => undef);
has program_version => (is => 'rw', init_arg => undef);
has command_line    => (is => 'rw', init_arg => undef);
has sort_order      => (is => 'rw', init_arg => undef);
has sam_version     => (is => 'rw', init_arg => undef);

# during constructor, read until first alignment and putback here
# want to use seek(), but doesn't work on stdin properly (?)
has putback         => (is => 'rw', init_arg => undef); 

# hash of sequence lengths. *Sequence names are uppercased*
has length          => (is => 'ro', default => sub { {} }, init_arg => undef); 

sub add_length{
    my ($self, $seqname, $length) = @_;
    croak "double adding $seqname into Sam::Parser? BUG, please report"
    if exists $self->length->{uc $seqname};
    $self->length->{uc $seqname} = $length;
}

# options
has convert_rc      => (is => 'ro', default => 0);
has skip_unmapped   => (is => 'ro', default => 1);


sub BUILD{
    my $self = shift;
    # read headers, putback first alignment line into $self->{putback}
    HEADER:
    while (defined(my $line = readline $self->filehandle)){
        chomp $line;
        if ($line =~ /^@/){
            $self->parse_header($line);
        }
        else{
            $self->putback($line);
            last HEADER;
        }
    }

    return $self;
}


# @HD	VN:1.0	SO:unsorted
# @SQ	SN:chr1	LN:30432563
# @SQ	SN:chr2	LN:19705359
# @SQ	SN:chr3	LN:23470805
# @SQ	SN:chr4	LN:18585042
# @SQ	SN:chr5	LN:26992728
# @SQ	SN:chrc	LN:154478
# @SQ	SN:chrm	LN:366924
# @PG	ID:Bowtie	VN:0.12.7	CL:"bowtie -S -f -B 1 -v 3 --best /wip/tools/genomes/AT/TAIR_reference.fas read.rc.fasta"
sub parse_header{
    my $self = shift;
    my $line = shift;

    my ($type_string, @parts) = split /\s+/, $line;

    # two letter code like @SQ
    my ($type) = $type_string =~ /^@(\w\w)/;

    die "can't parse header type" unless $type;

    my %header;

    # @SQ	[SN:chr5	LN:26992728] <= parse this part into %header
    for my $part (@parts) {
        if ($part =~ /(\w\w):(.*)/){
            my ($key, $value) = ($1, $2);
            $header{$key} = $value;
        }
    }

    if ($type eq 'HD' and defined($header{VN})){
        $self->sam_version($header{VN});
        if ($header{SO}){
            $self->sort_order($header{SO});
        }
    }
    elsif ($type eq 'SQ' and all { defined $header{ $_ } } qw/SN LN/){
        if ($self->convert_rc()){ 
            $header{SN} =~ s/^RC_//;
        }
        $self->add_length($header{SN}, $header{LN});
    }
    elsif ($type eq 'RG' and defined($header{ID})){
        # not sure what this is but leaving as stub
    }
    elsif ($type eq 'PG' and defined($header{ID})){
        $self->program_name($header{ID});
        if ($header{CL}){
            $self->command_line($header{CL});
        }
        if ($header{VN}){
            $self->program_version($header{VN});
        }
    }
    elsif ($type eq 'CO'){
        # ignore
    }
    else{
        die "can't parse header line";
    }
}

sub next{
    my $self = shift;

    # process putback
    if (defined (my $pb = $self->{putback})){
        $self->putback(undef);

        my $sam = Sam::samment->new($pb);
        if ($sam->mapped() || ! $self->skip_unmapped()){
            return $sam;
        }
    }

    while (defined(my $line = readline $self->file_handle)){
        my $sam = Sam::samment->new($line);
        if ($sam->mapped() || ! $self->skip_unmapped()){
            return $sam;
        }
    }

    return;
}

1;
