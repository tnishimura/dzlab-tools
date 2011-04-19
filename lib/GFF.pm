package GFF;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use Carp;
use autodie;    

has asterisk         => ( is => 'rw', isa => 'Bool', default => 0 );
# type constraints very slow? removed.
has sequence         => ( is => 'rw' );
has source           => ( is => 'rw' );
has feature          => ( is => 'rw' );
has start            => ( is => 'rw' );
has end              => ( is => 'rw' );
has score            => ( is => 'rw' );
has strand           => ( is => 'rw' );
has frame            => ( is => 'rw' );
has attribute_string => ( is => 'rw' );
has attr_hash => (
      traits    => ['Hash'],
      is        => 'ro',
      isa       => 'HashRef[Str]',
      handles   => {
          set_attribute  => 'set',
          get_attribute  => 'get',
          list_attribute => 'keys',
      },
      lazy_build => 1,
  );

sub _build_attr_hash{
    my $self = shift;
    my %accum;
    my $attrstr = $self->attribute_string;
    if (defined($attrstr) ){
        if ($attrstr=~s/\*$//){
            $self->asterisk(1);
            $self->attribute_string($attrstr);
        }
        for (split /;/, $attrstr){
            my ($key, $val) = split /=/, $_;
            $key =~ s/^\s+//;
            $key =~ s/\s+$//;

            if (defined $val){
                $val =~ s/^\s+//;
                $val =~ s/\s+$//;
                $accum{$key} = $val;
            }
            else {
                $accum{Note} = $key;
            }
        }
    }
    return \%accum;
}

sub to_string{
    my $self = shift;
    return 
    join "\t",
    map { ! defined $_ ? q{.} : $_ } 
    ($self->sequence, $self->source, $self->feature, $self->start, $self->end,
        $self->score,   $self->strand, $self->frame,   $self->attribute_string);
}

sub parse_locus{
    my ($self, $locus_tag) = @_;
    if ($self->get_attribute($locus_tag) =~ /([^.]+)\.([^.]+)/){
        return ($1,$2);
    } else{
        return;
    }
}

my %cols = map { $_ => 1 } qw/sequence source feature start end score strand frame attribute_string/;


sub get_column{
    my ($self, $colname) = @_;
    if (exists $cols{$colname}){
        return $self->$colname;
    } elsif ($colname =~ s/^\*//){
        my ($locus, $suffix) = $self->parse_locus($colname);
        return $locus;
    } else{
        return $self->get_attribute($colname);
    }
}

sub equals{
    my ($self, %against) = @_;
    my @columns = keys %against;
    foreach my $col (@columns) {
        my $x = $self->get_column($col);
        my $y = $against{$col};

        # both undefined, or both defined and equal
        if ((!defined $x && !defined $y) || (defined $x && defined $y && $x eq $y)){
            next;
        } else {
            return 0;
        }
    }
    return 1;
}

sub length{
    my $self = shift;
    return $self->end - $self->start + 1;
}

no Moose;
__PACKAGE__->meta->make_immutable;

#sub colname2num{ return $colmap{$_[0]} // croak "bad colmn name"; }

1;


=head1 NAME
 
GFF - GFF Moose Class
 
=head1 VERSION
 
This documentation refers to GFF version 0.0.1
 
=head1 SYNOPSIS
 
    use GFF
  
=head1 DESCRIPTION
 
<description>

=head1 SUBROUTINES/METHODS 

<exports>

=over

=item $gff->sequence

=item $gff->source  

=item $gff->feature 

=item $gff->start   

=item $gff->end     

=item $gff->score   

=item $gff->strand  

=item $gff->frame   

=item $gff->length

=item $gff->to_string

Format a GFF object for printing

=item $gff->get_column('sequence')

Like get_attribute, but also accepts qw/sequence source feature start end score
strand frame attribute_string/, as well as any attribute key.  If prefixed by a
'*', parse_locus called and the prefix is returned.  (So '*ID' on 'AT123.1'
returns 'AT123').

=cut

=item $gff->get_attribute('ID')

Get value of attribute in column 9.

=item $gff->set_attribute('ID') = 'AT12345'

Set value of attribute in column 9.

=item $gff->parse_locus('ID')

Returns the 'ID' attribute split by a dot ('.'), or undef.

=item $gff->equals({sequence => 'chr1', start => 1234})

Returns true if for every key $k in hash, $gff has $gff->get_column($k) eq $val.
For debugging mostly.

=back

=cut

