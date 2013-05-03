#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;
use IO::All;
use Pod::Usage;
use Getopt::Long;
use Parallel::ForkManager;
use LWP::Simple;
use XML::Simple;
use File::Path qw/make_path/;
use Path::Class;
use List::MoreUtils qw/all/;

my $result = GetOptions (
    "output-dir|d=s"    => \(my $output_dir = '.'),
    "num-processes|p=i" => \(my $num_processes = 4),
    "no-sra|n"          => \(my $no_sra),
    "dry"               => \(my $dry),
    "force|f"           => \(my $force),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || ! @ARGV || ! all { /^GSM\d+$/ } @ARGV);  

my $pm = Parallel::ForkManager->new($num_processes);

for my $gsmid (@ARGV) {
    my $basedir = "$output_dir/$gsmid";
    my @urls = get_all_supdata_urls_from_gsm($gsmid, $basedir) ;

    for my $url (@urls){
        $pm->start and next;
        my $basename = get_basename_from_url($url);
        my $file = "$basedir/$basename";

        if ($no_sra && $url =~ /\.sra$/i){
            say "[$gsmid] skipping b/c --no-sra $url";
            $pm->finish;
        }
        elsif ($dry){
            say "[$gsmid] skipping b/c --dry $url";
            $pm->finish;
        }
        elsif (-f $file and -s $file and ! $force){
            say "[$gsmid] skipping b/c it seems to already exist, and not doing --force $url";
            $pm->finish;
        }
        else{
            say "[$gsmid] grabbing $url";
            if (200 != getstore($url, $file)){
                say "couldn't grab $url?";
                exit 1;
            }
            $pm->finish; 
        }
    }
}
$pm->wait_all_children;

#######################################################################

# gsm is a geo sample. associated with it are supplemental files, which are
# either processed files, or raw files. raw files are SRX data sets, consisting
# of SRA files.

# parse out the part after the last "/"
sub get_basename_from_url{
    my $url = shift;
    return $url =~ s{^.*/([^/]+)}{$1}r;
}

# given id (gsm or srx), make the url, get meta file from ncbi, and return local name of meta file
sub get_meta{
    my ($id, $basedir) = @_;
    my $meta_file = "$basedir/$id.xml";
    my $meta_url = $id =~ /^SRX/ ? "http://www.ncbi.nlm.nih.gov/sra?term=$id&report=FullXml" 
                                 : "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$id&targ=self&form=xml&view=quick"; 
    say "grabbing $id meta file into $meta_file $meta_url";
    make_path($basedir);
    getstore($meta_url, $meta_file);
    return $meta_file;
}

sub get_all_sra_urls_from_srx{
    my ($srx, $basedir) = @_;
    my $meta_file = get_meta($srx, $basedir);
    my $meta_content = io($meta_file)->slurp;

    # <RUN alias="GSM543296_1" run_center="unspecified" accession="SRR052365" center_name="GEO" total_spots="15128663" total_bases="832076465" size="901689349" load_done="true" published="2010-06-18 11:51:24" is_public="true" cluster_name="public" static_data_available="1">
    my @all_srrids = $meta_content =~ m{accession="<span>(SRR\d{6})</span>"}xmg;
    my @urls;
    for my $srrid (@all_srrids) {
        # accession="SRR388657"
        # ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR388/SRR388657/SRR388657.sra
        if ($srrid =~ /(SRR\d\d\d)\d\d\d/){
            my $short = $1;
            push @urls, "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/$short/$srrid/$srrid.sra";
        }
        else{
            die "$srrid not a valid SRR number?";
        }
    }
    return @urls;
}

sub get_all_supdata_urls_from_gsm{
    my ($gsm, $basedir) = @_;
    my $meta_file = get_meta($gsm, $basedir);

    my $supdata = XMLin($meta_file, ForceArray => [qw/Supplementary-Data/])->{'Sample'}->{'Supplementary-Data'};
    # <Supplementary-Data type="SRA Experiment">
    # ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX111/SRX111004
    # </Supplementary-Data>

    my @urls;

    for my $d (@$supdata) {
        my $content = $d->{content}; 
        my $type    = $d->{type};
        $content =~ s/^\s+//; # strip off newlines
        $content =~ s/\s+$//;

        my $id = get_basename_from_url($content);

        if ($type  eq 'SRA Experiment'){
            # for SRX raw files, we need to parse the srx file for sra urls.
            push @urls, get_all_sra_urls_from_srx($id, $basedir);
        }
        else{
            # for processed files, the contents is just the url of the file
            push @urls, $content;
        }
    }
    return @urls;
}

=head1 geo-grab.pl 

Usage examples:

 geo-grab.pl [-d output-dir] [-p 4] [-f | --force] [--dry] GSM123456 GSM512341

=cut


__DATA__
sample gsm meta file (relevant parts):

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM847326&targ=self&form=xml&view=quick

<?xml version="1.0" encoding="UTF-8" standalone="no"?>

<MINiML
   xmlns="http://www.ncbi.nlm.nih.gov/geo/info/MINiML"
   xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
   xsi:schemaLocation="http://www.ncbi.nlm.nih.gov/geo/info/MINiML http://www.ncbi.nlm.nih.gov/geo/info/MINiML.xsd"
   version="0.5.0" >

  <Sample iid="GSM847326">
    <Contact-Ref ref="contrib1" />
    <Supplementary-Data type="TXT">
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM847nnn/GSM847326/suppl/GSM847326_wdleaf_lan1_bowtie.txt.gz
    </Supplementary-Data>
    <Supplementary-Data type="SRA Experiment">
    ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX111/SRX111004
    </Supplementary-Data>
    <Relation type="SRA" target="http://www.ncbi.nlm.nih.gov/sra?term=SRX111004" />
  </Sample>

</MINiML>

