package Sam::Parser;
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Hash::Util qw/lock_keys unlock_keys/;
use parent 'ParserNoMoose';
use List::MoreUtils qw/all/;

sub new {
    my $class = shift;
    my %opt = @_;

    my $convert_rc = delete $opt{convert_rc};
    my $skip_unmapped = delete $opt{skip_unmapped};

    my $self = $class->SUPER::new(%opt);

    unlock_keys(%$self);

    $self->{length}          = {};
    $self->{rc_sequence}     = {};    # rc_sequence caches the names of the ^RC_ chromosomes
    $self->{program_name}    = undef;
    $self->{program_version} = undef;
    $self->{command_line}    = undef;
    $self->{sort_order}      = undef;
    $self->{sam_version}     = undef;
    $self->{convert_rc}      = $convert_rc // 0;
    $self->{skip_unmapped}   = $skip_unmapped // 1;
    $self->{putback}         = undef; # during constructor, read until first alignment and putback here
                                      # want to use seek(), but doesn't work on stdin properly (?)

    lock_keys(%$self);

    # read headers, putback first alignment line into $self->{putback}
    HEADER:
    while (defined(my $line = readline $self->{handle})){
        chomp $line;
        if ($line =~ /^@/){
            $self->parse_header($line);
        }
        else{
            $self->{putback} = $line;
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

    my ($type) = $type_string =~ /^@(\w\w)/;

    die "can't parse header type" unless $type;

    my %header;

    for my $part (@parts) {
        if ($part =~ /(\w\w):(.*)/){
            my ($key, $value) = ($1, $2);
            $header{$key} = $value;
        }
    }

    if ($type eq 'HD' and defined($header{VN})){
        $self->{sam_version} = $header{VN};
        if ($header{SO}){
            $self->{sort_order} = $header{SO};
        }
    }
    elsif ($type eq 'SQ' and all { defined $header{ $_ } } qw/SN LN/){
        $self->{length}{$header{SN}} = $header{LN};
        if ($self->{convert_rc} and $header{SN} =~ /^RC_/){
            $self->{rc_sequence}{$header{SN}} = 1;
        }
    }
    elsif ($type eq 'RG' and defined($header{ID})){
        # not sure what this is but leaving as stub
    }
    elsif ($type eq 'PG' and defined($header{ID})){
        $self->{program_name} = $header{ID};
        if ($header{CL}){
            $self->{command_line} = $header{CL};
        }
        if ($header{VN}){
            $self->{program_version} = $header{VN};
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
        undef $self->{putback};
        my $align = $self->parse_alignment($pb);
        if ($align->{mapped} || ! $self->{skip_unmapped}){
            return $align;
        }
    }

    while (defined(my $line = readline $self->{handle})){
        chomp $line;
        my $align = $self->parse_alignment($line);
        if ($align->{mapped} || ! $self->{skip_unmapped}){
            return $align;
        }
    }

    return;
}

1;
