package Sam::Parser;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Hash::Util qw/lock_keys unlock_keys/;
use List::MoreUtils qw/all/;
use Moose;
use Sam::Alignment;

extends 'Parser';

has '+filename_or_handle' => ( required => 0,);

# from headers
has program_name    => (is => 'rw', init_arg => undef);
has program_version => (is => 'rw', init_arg => undef);
has command_line    => (is => 'rw', init_arg => undef);
has sort_order      => (is => 'rw', init_arg => undef);
has sam_version     => (is => 'rw', init_arg => undef);

# during constructor, read until first alignment and putback here
# want to use seek(), but doesn't work on stdin properly (?)
has putback => (is => 'rw', init_arg => undef); 

# hash of sequence lengths. *Sequence names are uppercased*
has length  => (is => 'ro', default => sub { {} }, init_arg => undef); 

sub add_length{
    my ($self, $seqname, $length) = @_;

    if (exists $self->length->{uc $seqname} and $self->length->{uc $seqname} != $length){
        my $incumbent = $self->length->{uc $seqname};
        croak "double adding conflicting lengths for $seqname into Sam::Parser? New length = $length, incumbent = $incumbent. BUG, please report"
    }
    $self->length->{uc $seqname} = $length;
}

has header_lines => (
    is      => 'ro',
    default => sub { [] },
);

sub header{
    my $self = shift;

    my @header_lines = @{$self->header_lines};
    chomp @header_lines;
    return join "\n", @header_lines;
}

# options
has convert_rc      => (is => 'ro', default => 0);
has skip_unmapped   => (is => 'ro', default => 1);

sub BUILD{
    my $self = shift;
    return if (! defined $self->filename_or_handle);

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

# return the next valid Sam::Alignment object, or undef
sub next{
    my $self = shift;

    # process putback
    if (defined (my $pb = $self->putback())){
        $self->putback(undef);

        my $sam = Sam::Alignment->new($pb, $self->length(), $self->convert_rc);
        if ($sam->mapped() || ! $self->skip_unmapped()){
            return $sam;
        }
    }

    return if eof($self->filehandle);

    while (defined(my $line = readline $self->filehandle)){
        my $sam = Sam::Alignment->new($line, $self->length(), $self->convert_rc);
        if ($sam->mapped() || ! $self->skip_unmapped()){
            return $sam;
        }
    }

    return;
}

# for use when filename_or_handle is omitted and you want to use as push parser
# returns: Sam::Alignment object OR header line OR undef (when given header,
# unmapped line when skipping unmapped)
sub push{
    my ($self, $line) = @_;
    if ($line =~ /^@/){
        $self->parse_header($line);
        return $line;
    }
    else{
        my $sam = Sam::Alignment->new($line, $self->length(), $self->convert_rc);
        if ($sam->mapped() || ! $self->skip_unmapped()){
            return $sam;
        }
    }
    return;
}

# given a header line, fills in Sam::Parser's parameters according and 
# appends line to header_lines
# sample:
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

    CORE::push @{$self->header_lines}, $line;

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
            if ($header{SN} =~ s/^RC_//){
                # don't want to have RC_chr* in header if we're converting
                pop @{$self->header_lines};
            }
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

1;
