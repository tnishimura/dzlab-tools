#!/usr/bin/env perl
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use IO::File;
use IO::Handle;
use XML::Simple;
use LWP::Simple;
use Pod::Usage;
use Getopt::Long;
use Path::Class;
use LWP::Simple;
use File::Path qw/make_path/;


my $result = GetOptions (
    "output-dir|d=s" => \(my $output_dir),
    "gsm=s"          => \(my $gsm),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || !$output_dir || ! $gsm);  

$output_dir = dir($output_dir);
make_path($output_dir);

my $miniml = $output_dir->file("$gsm.xml")->stringify;
my $miniml_url = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$gsm&targ=self&form=xml&view=quick";

warn "grabbing $miniml_url into $output_dir";
getstore($miniml_url, $miniml);

my $supdata = XMLin($miniml, ForceArray => [qw/Supplementary-Data/])->{'Sample'}->{'Supplementary-Data'};

for my $d (@$supdata) {
    my $content = $d->{content};
    my $type    = $d->{type};
    $content =~ s/^\s+//;
    $content =~ s/\s+$//;

    my $basename = $content =~ s{^.*/([^/]+)}{$1}r;
    my $file = $type eq 'SRA Experiment' ? $output_dir->file("$basename.sra") : $output_dir->file("$basename");

    warn "grabbing $content into $file";
    getstore($content, $file->stringify);
}

=head1 geo-grab.pl 

Usage examples:

 geo-grab.pl -d output-dir -gsm GSM123456

=cut


# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM847326&targ=self&form=xml&view=quick
# 
# <?xml version="1.0" encoding="UTF-8" standalone="no"?>
# 
# <MINiML
#    xmlns="http://www.ncbi.nlm.nih.gov/geo/info/MINiML"
#    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
#    xsi:schemaLocation="http://www.ncbi.nlm.nih.gov/geo/info/MINiML http://www.ncbi.nlm.nih.gov/geo/info/MINiML.xsd"
#    version="0.5.0" >
# 
#   <Sample iid="GSM847326">
#     <Contact-Ref ref="contrib1" />
#     <Supplementary-Data type="TXT">
# ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM847nnn/GSM847326/suppl/GSM847326_wdleaf_lan1_bowtie.txt.gz
#     </Supplementary-Data>
#     <Supplementary-Data type="SRA Experiment">
#     ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX111/SRX111004
#     </Supplementary-Data>
#     <Relation type="SRA" target="http://www.ncbi.nlm.nih.gov/sra?term=SRX111004" />
#   </Sample>
# 
# </MINiML>
