#!/usr/bin/env perl
package main;
use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config gnu_getopt);
use Pod::Usage;
use File::Basename;
use PPI;
use List::Util 'max';
use Algorithm::Permute qw(permute);
use Statistics::Basic qw(:all nofill);

my $INH  = *ARGV;
my $ERRH = *STDERR;
my $OUTH = *STDOUT;
GetOptions(
    \%ARGV,
    'input|i=s',
    'output|o=s',
    'debug:i',
    'quiet'   => sub { $ARGV{quiet}   = 1; no diagnostics;  no warnings },
    'verbose' => sub { $ARGV{verbose} = 1; use diagnostics; use warnings },
    'version' => sub {
        pod2usage
          -sections => [ 'VERSION', 'REVISION' ],
          -verbose  => 99;
    },
    'license' => sub {
        pod2usage
          -sections => [ 'AUTHOR', 'COPYRIGHT' ],
          -verbose  => 99;
    },
    'usage' => sub {
        pod2usage
          -sections => ['SYNOPSIS'],
          -verbose  => 99;
    },
    'help'   => sub { pod2usage -verbose => 1 },
    'manual' => sub { pod2usage -verbose => 2 },
) or pod2usage( -verbose => 1 );

IO:
{

    # We use the ARGV magical handle to read input
    # If user explicitly sets -i, put the argument in @ARGV
    if ( exists $ARGV{input} ) {
        unshift @ARGV, $ARGV{input};
    }

    # Allow in-situ arguments (equal input and output filenames)
    # FIXME: infinite loop. Why?
    if (    exists $ARGV{input}
        and exists $ARGV{output}
        and $ARGV{input} eq $ARGV{output} )
    {
        croak "Bug: don't use in-situ editing (same input and output files";
        open $INH, q{<}, $ARGV{input}
          or croak "Can't read $ARGV{input}: $!";
        unlink $ARGV{input};
    }

    # Redirect STDOUT to a file if so specified
    if ( exists $ARGV{output} and q{-} ne $ARGV{output} ) {
        open $OUTH, q{>}, $ARGV{output}
          or croak "Can't write $ARGV{output}: $!";
    }
}

my %data;
foreach my $file (@ARGV) {

    # Load a document from a file
    my $document = PPI::Document->new($file);

    # Strip out comments
    $document->prune('PPI::Token::Comment');

    # $document->normalized;
    # $document->index_locations;

    # Find all the named subroutines
    my $sub_nodes =
	$document->find( sub { $_[1]->isa('PPI::Statement::Sub') and $_[1]->name });

    next unless $sub_nodes;
    my %sub_names = map {
        $_->name => {
            'code'   => $_->content,
	    'tokens' => token_frequency($_)
	}
    } @$sub_nodes;

    my $file_name = fileparse $file;
    $data{$file_name} = \%sub_names;
}




#Set up @file_names and @sub_names for easy reference
my @file_names = sort keys %data;
my @sub_names;
foreach my $sub_hash ( values %data ) {
    foreach my $sub_name ( keys %{$sub_hash} ) {
        push @sub_names, $sub_name;
    }
}
#remove duplicates and sort
my %hash = map { $_, 1 } @sub_names;
@sub_names = sort keys %hash;



#Put a third key into each subroutine hash, pointing back to the hash for each file in common.
#Each element in the array is a hash.
foreach my $sub_name (@sub_names) {
    my @file_names_by_sub;
    foreach my $file_name (@file_names) {
        if ( exists $data{$file_name}->{$sub_name} ) {
            push @file_names_by_sub, $file_name;
        }
    }
    foreach my $file_name_by_sub (@file_names_by_sub) {
        foreach my $file_name (@file_names) {
            $data{$file_name}->{$sub_name}->{'other files'}->{$file_name_by_sub} = $data{$file_name} 
	       if exists $data{$file_name}->{$sub_name}
	          and $file_name ne $file_name_by_sub;   #keeps it from pointing to itself, and from creating new keys.

	}
    }
}

{
    local $Data::Dumper::Maxdepth = 4;
#    print Dumper keys %{$data{'parse_bowtie.pl'}{'index_fasta'}{'other files'}};
 #   exit;
}


##should I unify %data and %cc_by_sub_name? 



#cross correlation part 
#do it for different names, too.
my %cc_by_sub_name;

foreach my $sub_name (@sub_names) {
    my %code_of_sub_variations;
    foreach my $file ( keys %data ) {
        foreach my $my_sub_name ( keys %{ $data{$file} } ) {
            if ( $sub_name eq $my_sub_name) {
                $code_of_sub_variations{$file} =
                  $data{$file}->{$my_sub_name}->{'tokens'};
            }
        }
    }
    my $sub_cc = cc_all_combos(%code_of_sub_variations);
    $cc_by_sub_name{$sub_name} = $sub_cc
	if %$sub_cc;    #gets rid of subs used in only one file.
}

#die Dumper \%cc_by_sub_name;


auto_module(\%cc_by_sub_name, "/home/jgraff/workspace/bisulfite/trunk/module_test.pm");

#my %thing = subs_of_interest(1, \%cc_by_sub_name, @sub_names);
#my %thing = subs_of_interest(@sub_names,);


#print Dumper \%thing; exit;
#print Dumper \%cc_by_sub_name;



##More stuff to do: add mode switches for table, subs of interest, and implement auto-generation of modules
#subs of interest thing should take a tolerance (0 - 1) and return all files/subs that exceed that number. Output: lists of subnames/filenames.
#auto-gen modules: find cases where correlation is perfect, and put those into a single module.
#Turn the whole thing into a module, make methods instead of modes.



#I should examine the output more closely and make sure it matches up with the numbers.
sub auto_module {
    my ($cc_by_sub_name, $module_file) = @_;
    my %cc_by_sub_name = %$cc_by_sub_name;
    my %cced_sub_code = _get_perfect_cc_subs(%cc_by_sub_name);	    	       	
    unlink $module_file;
    _create_module($module_file, sort keys %cced_sub_code)
	unless -e $module_file;
    _add_to_module($_, $module_file) foreach sort values %cced_sub_code;	          
}
 

sub _get_perfect_cc_subs {
    my (%cc_by_sub_name) = @_;
    my %cced_sub_code;
    for my $sub_name (keys %cc_by_sub_name) {
	my @files = keys %{ $cc_by_sub_name{$sub_name} };
	my $file1 = $files[0];
	for my $cc (values %{ $cc_by_sub_name{$sub_name}{$file1} }) {
	    if ($cc == 1) {
		my $code = $data{$file1}{$sub_name}{'code'};
		$cced_sub_code{$sub_name} = $code;
		last;
	    }
	}
    }
    return %cced_sub_code;
}



sub _add_to_module {
    my ($code, $module) = @_;
    open my $MODULE, '>>', $module or croak "Can't open $module: $!";
    print $MODULE "\n$code\n";
    close $MODULE or  croak "Can't close $module: $!";
}

sub _create_module {
    my ($file, @sub_names) = @_;
    open my $FILE, '>', $file or croak "Can't open $file: $!";
    print $FILE 
	"package Utils; \n",
	"use Export;\n",
	"our \@EXPORT_OK = qw(@sub_names);\n"; #basic module stuff 
    
    close $FILE or croak "Can't close $file: $!";
}


#####This needs improvement! The groups are repeated w/ diff permutations. Possibly rethink approach.

#returns a hash, keys are sub names, values is an array of files.
#I can easily change this back to an array.
sub subs_of_interest {
    my ($tolerance, $cc_by_sub_name, @sub_names) = @_;
    my %cc_by_sub_name = %$cc_by_sub_name;
    my %subs_of_interest;
    for my $sub_name (keys %cc_by_sub_name) {		
	my @groups;
	my @done;
	for my $file1 (keys %{$cc_by_sub_name{$sub_name}}) {	    	
	    my @files = ($file1);	    
	    for my $file2 (keys %{$cc_by_sub_name{$sub_name}{$file1}}) {        
		my $cc = $cc_by_sub_name{$sub_name}{$file1}{$file2};
		push @files, $file2
		    unless $cc < $tolerance or grep /$file2/, @done;   		    
		push @done, $file1;
	    }
	    push @groups, \@files
		if scalar @files > 1; #if there's only one element, there are no matches for the given tolerance.
	}
	$subs_of_interest{$sub_name} = \@groups;
	
    }
    return %subs_of_interest;
}





# Takes a hash representing a single subroutine. keys are files, values are hashes of the token frequencies.
# Returns a hash with keys being a filename pointing to a hash whose keys are another filename and whose value is the token comparison for the sub between those two files.
# So to retieve the difference between file2 and file3, do $hash{file2}{file3} or $hash{file3}{file2}.
sub cc_all_combos {
    my (%hash) = @_;
    my @files = keys %hash;
    my %ret_hash;
    my @combinations;

    while (my $file1 = shift @files) {
	foreach my $file2 (@files) {
	    my @combo = ($file1, $file2);
	    push @combinations, \@combo;
	}
    }

    for my $combo (@combinations) {
	my $file1 = @{$combo}[0];
	my $file2 = @{$combo}[1];
	$ret_hash{$file1}->{$file2} = cross_correlation($hash{$file1}, $hash{$file2});
	$ret_hash{$file2}->{$file1} = cross_correlation($hash{$file1}, $hash{$file2});
    }
    return \%ret_hash;
}



sub normalized_hd {
  my ($k, $l) = @_;
  my $xor = $k ^ $l;
  my $hamming_distance = $xor =~ tr/\0//c;
  my $max_length = max length $k, length $l;
  return sprintf ("%g", 1 - $hamming_distance / $max_length)
      unless $max_length == 0;
  return -1;
}



sub cross_correlation {
    my ($tokens1, $tokens2) = @_; 
    
    my %all_tokens = map { $_ => 0 } (keys %$tokens1, keys %$tokens2);
    my (@vec1, @vec2);

    for my $token (keys %all_tokens) { 

	$tokens1->{$token} = 0 unless $tokens1->{$token};
	$tokens2->{$token} = 0 unless $tokens2->{$token};

	push @vec1, $tokens1->{$token};
	push @vec2, $tokens2->{$token};
    }	

    my $cc = correlation(vector(@vec1), vector(@vec2));
    return "$cc";
}



#takes a sub_node object, returns a ref to a hash of token frequencies
sub token_frequency {
    my ($sub_node) = @_;
    
    my @tokens = $sub_node->find( 
	sub {
	    $_[1]->isa('PPI::Token')
		and not $_[1]->isa('PPI::Token::Whitespace')
	} );
    my %freqs;
    $freqs{ join "\t", ref ($_), $_->content }++ for @{$tokens[0]}; 
    return \%freqs;    
}



sub norm {
    my ( $string, @strings ) = @_;
    my $longest = max( map { length } @strings );
    until ( length $string == $longest ) {
        $string = "$string ";
    }
    return $string;
}


sub print_table {
##print out a table showing which files contain which subroutines.
    my (@sub_names, @file_names, %data) = @_;
    print join( "\t", q{}, @sub_names ), "\n";
    
    foreach my $file_name (@file_names) {
	print $file_name, "\t";
	foreach my $sub_name (@sub_names) {
	    if ( $data{$file_name}->{$sub_name} ) {
		print q{*};
	    }
	    else { print q{}; }
	    print "\t";
	}
	print "\n";
    }
}



__DATA__

__END__
=head1 NAME

 APerlyName.pl - Short description

=head1 SYNOPSIS

 APerlyName.pl [OPTION]... [FILE]...

=head1 DESCRIPTION

 Long description

=head1 OPTIONS

 -i, --input       filename to read from                            (STDIN)
 -o, --output      filename to write to                             (STDOUT)
     --debug       print additional information
     --verbose     print diagnostic and warning messages
     --quiet       print no diagnostic or warning messages
     --version     print current version
     --license     print author's contact and copyright information
     --help        print this information
     --manual      print the plain old documentation page

=head1 VERSION

 0.0.1

=head1 REVISION

 $Rev: $:
 $Author: $:
 $Date: $:
 $HeadURL: $:
 $Id: $:

=head1 AUTHOR

 Pedro Silva <pedros@berkeley.edu/>
 Zilberman Lab <http://dzlab.pmb.berkeley.edu/>
 Plant and Microbial Biology Department
 College of Natural Resources
 University of California, Berkeley

=head1 COPYRIGHT

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
