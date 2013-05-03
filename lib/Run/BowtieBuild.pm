package Run::BowtieBuild;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use Carp;
use autodie;
use List::MoreUtils qw/notall/;
use Params::Validate qw/:all/;

use FastaUtil;
use Launch qw/cast/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(bowtie_build);

=head2 

 my ($converted_file_name, $index_prefix_name, @ebwt_files) = bowtie_build(
     file => 'reference.fas', 
     noref => BOOL, 
     bs => ['c2t' | 'g2a'], 
     rc => BOOL, 
     force => BOOL,
     colorspace => BOOL,
 )

Returns ($converted_file_name, $index_prefix_name, @ebwt_files);

=cut
sub bowtie_build{
    my %opt = validate(@_, {
            file => {
                type => SCALAR,
                callbacks => {
                    'file exists' => sub { -f shift },
                },
            }, 
            noref => { default => 0, },
            index => { optional => 1, },
            bs => {
                type => SCALAR,
                callbacks => {
                    'bs should be c2t or g2a' => sub {
                        my $bstype = shift;
                        $bstype eq 'c2t' || $bstype eq 'g2a'
                    }
                },
                optional => 1,
            },
            rc => {
                default => 1, 
                optional => 1
            },
            force => 0,
            version => {
                optional => 1,
                default => 1, 
                regex => qr/^1|2$/,
            },
            colorspace => 0,
        });

    my $executable = $opt{version} == 1 ? 'bowtie-build' : 'bowtie2-build';
    my $suffix     = $opt{version} == 1 ? 'ebwt' : 'bt2';

    if ($opt{bs} && $opt{rc}){
        my $bsrc_file = "$opt{file}.$opt{bs}";
        if ($opt{force} || ! -f $bsrc_file){
            bsrc_fasta_on_disk($opt{bs}, $opt{file}, $bsrc_file);
        }
        $opt{file} = $bsrc_file;
    }
    elsif ($opt{bs} && ! $opt{rc}){ # naming conflict if you want to do both....
        my $bs_file = "$opt{file}.$opt{bs}";
        if (! -f $bs_file){
            bs_fasta_on_disk($opt{bs}, $opt{file}, $bs_file);
        }
        $opt{file} = $bs_file;
    }

    if (! $opt{index}){
        if ($opt{colorspace}){
            $opt{index} = $opt{file} . ".cs";
        }
        else{
            $opt{index} = $opt{file};
        }
    }

    my @expected_files = map { "$opt{index}.$_" } ("1.$suffix","2.$suffix","3.$suffix","4.$suffix", $opt{noref} ? () : ("3.$suffix", "4.$suffix"));

    my $norefarg = $opt{noref} ? '--noref' : '';
    if ($opt{force} or notall { -f && -s } @expected_files){
        cast "$executable $norefarg $opt{file} $opt{index}";
    }
    return ($opt{file}, $opt{index}, @expected_files);
}

1;

