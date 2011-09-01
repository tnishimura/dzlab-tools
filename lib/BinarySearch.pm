package BinarySearch;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use Scalar::Util qw/looks_like_number/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(greatest_lower least_upper);
our @EXPORT = qw();

sub _extended_cmp{
    my ($x, $y) = @_;
    my ($lcx, $lcy) = (lc $x, lc $y);
    my ($x_inf, $y_inf) = ($lcx eq 'inf' || $lcx eq '+inf', $lcy eq 'inf' || $lcy eq '+inf');
    my ($x_ninf, $y_ninf) = ($lcx eq '-inf', $lcy eq '-inf');
    # b/c looks_like_number('inf') is true...
    my ($x_lln, $y_lln) = (!$x_inf && !$x_ninf && looks_like_number($x), !$y_inf && !$y_ninf && looks_like_number($y));

    if ($x_lln && $y_lln){
        return $x <=> $y;
    }
    elsif ($x_lln){
        if ($y_inf){ return -1; }
        elsif ($y_ninf){ return 1; }
    }
    elsif ($y_lln){
        if ($x_inf) { return 1; }
        elsif ($x_ninf){ return -1; }
    }
    elsif ($x_ninf && $y_inf){ return -1; }
    elsif ($x_inf && $y_ninf){ return 1; }
    else {
        return;
    }
}

sub _extended_lt{ return _extended_cmp($_[0], $_[1]) == -1; }
sub _extended_gt{ return _extended_cmp($_[0], $_[1]) == 1; }
sub _extended_le{
    my $res = _extended_cmp($_[0], $_[1]);
    return $res == -1 || $res == 0;
}
sub _extended_ge{
    my $res = _extended_cmp($_[0], $_[1]);
    return $res == 1 || $res == 0;
}

sub _greatest_lower_linear{
    my ($array, $target) = @_;
    my @list = ('-inf', @$array, 'inf');
    for my $index (0 .. $#list - 1) {
        if (_extended_lt($list[$index],$target) && _extended_le($target, $list[$index+1])){
            my $val = $list[$index] ;
            if ($val eq '-inf'){
                return;
            }
            else{
                return $val;
            }
        }
    }
}
sub _least_upper_linear{
    my ($array, $target) = @_;
    my @list = ('-inf', @$array, 'inf');
    for my $index (0 .. $#list - 1) {
        if (_extended_le($list[$index],$target) && _extended_lt($target, $list[$index+1])){
            my $val = $list[$index+1];
            if ($val eq 'inf'){
                return;
            }
            else {
                return $val;
            }
        }
    }
}

sub greatest_lower{
    my ($array, $target) = @_;
    my @list = ('-inf', @$array, 'inf');
    my $result_index = _gl_helper(\@list, $target, 0, $#list-1);

    if ($result_index == 0){
        return;
    }
    return $list[$result_index];
}

sub _gl_helper{
    my ($array, $target, $low, $high) = @_;

    my $mid = $low + int(($high - $low)/2);
    my $mid_val = $array->[$mid];
    my $forward_val = $array->[$mid + 1];
    if (_extended_lt($mid_val, $target) && _extended_le($target, $forward_val)){
        return $mid;
    }
    elsif (_extended_ge($mid_val, $target)){
        return _gl_helper($array, $target, $low, $mid-1);
    }
    elsif (_extended_lt($forward_val, $target)){
        return _gl_helper($array, $target, $mid+1, $high);
    }
    else {
        croak "shouldn't be here";
    }
}

sub least_upper{
    my ($array, $target) = @_;
    my @list = ('-inf', @$array, 'inf');
    my $result_index = _lu_helper(\@list, $target, 0, $#list -1);

    if ($result_index == $#list){
        return;
    }
    return $list[$result_index];
}

sub _lu_helper{
    my ($array, $target, $low, $high) = @_;

    my $mid = $low + int(($high - $low)/2);
    my $mid_val = $array->[$mid];
    my $back_val = $array->[$mid - 1];
    if (_extended_le($back_val, $target) && _extended_lt($target, $mid_val)){
        return $mid;
    }
    elsif (_extended_le($mid_val, $target)){
        return _lu_helper($array, $target, $mid+1, $high);
    }
    elsif (_extended_lt($target, $back_val)){
        return _lu_helper($array, $target, $low, $mid - 1);
    }
    else {
        die "shouldn't be here";
    }

}
1;

