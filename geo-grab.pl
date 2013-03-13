#!/usr/bin/env perl
package GEO{
    use strict;
    use warnings;
    use 5.010_000;
    use Mouse;
    use Carp;
    use LWP::Simple;
    use IO::All;
    use File::Path qw/make_path/;

    sub BUILD{
        my $self = shift;
        make_path($self->base_directory);
        getstore($self->meta_url, $self->meta_file);
    }

    has 'id'        => (required => 1, is => 'ro');
    has 'base_directory' => (required => 1, is => 'ro');

    sub meta_url{
        my $self = shift;
        my $id = $self->id;

        if ($id =~ /^SRX/){
            return "http://www.ncbi.nlm.nih.gov/sra?term=$id&report=FullXml";
        }
        else{
            return "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$id&targ=self&form=xml&view=quick";
        }
    }

    sub meta_file{
        my $self = shift;
        my $dir = $self->base_directory();
        my $id = $self->id();
        return "$dir/$id.xml";
    }

    sub meta_content{
        my $self = shift;
        return io($self->meta_file)->slurp;
    }
}
package GEO::SRX{
    use strict;
    use warnings;
    use 5.010_000;
    use Mouse;
    use Carp;
    use XML::Simple;
    use LWP::Simple;
    use IO::All;

    extends 'GEO';

    # http://www.ncbi.nlm.nih.gov/sra?term=SRX021424&report=FullXml
    # return list of url
    sub all_sra{
        my $self = shift;
        my @srrid = $self->meta_content =~ m{accession="<span>(SRR\d{6})</span>"}xmg;
        return map { 
            make_srr_url($_) 
        } @srrid;
        # <RUN alias="GSM543296_1" run_center="unspecified" accession="SRR052365" center_name="GEO" total_spots="15128663" total_bases="832076465" size="901689349" load_done="true" published="2010-06-18 11:51:24" is_public="true" cluster_name="public" static_data_available="1">
    }

    sub make_srr_url{
        # accession="SRR388657"
        # ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR388/SRR388657/SRR388657.sra
        my $srrid = shift;
        if ($srrid =~ /(SRR\d\d\d)\d\d\d/){
            my $short = $1;
            return "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/$short/$srrid/$srrid.sra";
        }
        else{
            die "$srrid not a valid SRR number?";
        }
        # ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR388/SRR388657/SRR388657.sra
    }
}
package GEO::GSM{
    use strict;
    use warnings;
    use 5.010_000;
    use Mouse;
    use Carp;
    use XML::Simple;
    use LWP::Simple;
    use File::Path qw/make_path/;
    use Path::Class;

    extends 'GEO';

    sub gsm_directory{
        my $self = shift;
        my $subdir = dir($self->base_directory())->subdir($self->id)->stringify;
        make_path($subdir);
        return $subdir;
    }

    sub _supdata{
        my $self = shift;
        my $supdata = XMLin($self->meta_file, ForceArray => [qw/Supplementary-Data/])->{'Sample'}->{'Supplementary-Data'};
        # <Supplementary-Data type="SRA Experiment">
        # ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX111/SRX111004
        # </Supplementary-Data>

        for my $d (@$supdata) {
            my $content = $d->{content}; 
            # my $type    = $d->{type};
            $d->{content} =~ s/^\s+//;
            $d->{content} =~ s/\s+$//;
        }
        return @$supdata
    }

    sub all_srx{
        my $self = shift;
        my @supdata = $self->_supdata;

        my @accum;
        for my $d (@supdata) {
            my $content = $d->{content}; 
            my $type    = $d->{type};
            my $basename = $content =~ s{^.*/([^/]+)}{$1}r;

            if ($type  eq 'SRA Experiment'){
                push @accum, GEO::SRX->new(id => $basename, base_directory => $self->base_directory);
            }
        }
        return @accum;
    }

    # grab all sra's from all srx's.
    sub all_sra{
        my $self = shift;
        my @srx = $self->all_srx;
        map 
        {
            $_->all_sra()
        } @srx;
    }

    sub all_processed{
        my $self = shift;
        my @supdata = $self->_supdata;
        my @accum;
        for my $d (@supdata) {
            my $content = $d->{content}; 
            my $type    = $d->{type};
            my $basename = $content =~ s{^.*/([^/]+)}{$1}r;

            if ($type ne 'SRA Experiment'){
                push @accum, $content; 
            }
        }
        return @accum;
    }

    sub all_supplementary_files{
        my $self = shift;
        return sort $self->all_sra(), $self->all_processed();
    }
}


package GEO::GSE{
    use strict;
    use warnings;
    use 5.010_000;
    use Mouse;
    use Carp;
    use XML::Simple;
    use LWP::Simple;

    extends 'GEO';

    sub all_gsm{
        my $self = shift;
        $self->meta_file;

        my $samples = XMLin($self->meta_file, ForceArray => [qw/Sample/])->{Sample};
        return map 
        {
            my $gsmid = $_->{iid};
            GEO::GSM->new(id => $gsmid, base_directory => $self->base_directory() . "/$gsmid");
        } @$samples;
    }
}


package main;
use strict;
use warnings FATAL => "all";
use 5.010_000;
use Data::Dumper;
use autodie;
use Pod::Usage;
use Getopt::Long;
use Parallel::ForkManager;
use LWP::Simple;

my $result = GetOptions (
    "output-dir|d=s" => \(my $output_dir = '.'),
    "num-processes|p=i" => \(my $num_processes = 0),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (!$result || !$output_dir || @ARGV==0);  

my $pm = Parallel::ForkManager->new($num_processes);

my @gsm = @ARGV;

for my $id (@gsm) {
    say "[$id] starting";
    my $gsmdir = "$output_dir/$id/";
    my $gsm = GEO::GSM->new(id => $id, base_directory => $gsmdir);
    for my $url ($gsm->all_supplementary_files) {
        $pm->start and next;
        my $basename = $url =~ s{^.*/([^/]+)}{$1}r;

        say "[$id] grabbing $url";
        my $rc = getstore($url, "$gsmdir/$basename");

        if ($rc != 200){
            die "couldn't grab $url?";
        }
        $pm->finish; 
    }
}
$pm->wait_all_children;

=head1 geo-grab.pl 

Usage examples:

 geo-grab.pl -d output-dir [-p 4] GSE76543 GSM123456

=cut


# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34318&targ=self&form=xml&view=quick
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

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34318&targ=self&form=xml&view=quick
# <?xml version="1.0" encoding="UTF-8" standalone="no"?>
# <MINiML
#    xmlns="http://www.ncbi.nlm.nih.gov/geo/info/MINiML"
#    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
#    xsi:schemaLocation="http://www.ncbi.nlm.nih.gov/geo/info/MINiML http://www.ncbi.nlm.nih.gov/geo/info/MINiML.xsd"
#    version="0.5.0" >
# 
#   <Contributor iid="contrib1">
#     <Person><First>Yufeng</First><Last>Wu</Last></Person>
#     <Email>ywu38@wisc.edu</Email>
#     <Laboratory>Jiming Jiang</Laboratory>
#     <Organization>UW-Madison</Organization>
#     <Address>
#       <Line>1575 Linden Drive</Line>
#       <City>Madison</City>
#       <State>WI</State>
#       <Zip-Code>53706</Zip-Code>
#       <Country>USA</Country>
#     </Address>
#   </Contributor>
# 
#   <Contributor iid="contrib2">
#     <Person><First>Wenli</First><Last>Zhang</Last></Person>
#   </Contributor>
# 
#   <Contributor iid="contrib3">
#     <Person><First>Tao</First><Last>Zhang</Last></Person>
#   </Contributor>
# 
#   <Contributor iid="contrib4">
#     <Person><First>Jiming</First><Last>Jiang</Last></Person>
#   </Contributor>
# 
#   <Database iid="GEO">
#     <Name>Gene Expression Omnibus (GEO)</Name>
#     <Public-ID>GEO</Public-ID>
#     <Organization>NCBI NLM NIH</Organization>
#     <Web-Link>http://www.ncbi.nlm.nih.gov/geo</Web-Link>
#     <Email>geo@ncbi.nlm.nih.gov</Email>
#   </Database>
# 
#   <Sample iid="GSM847326">
#     <Accession database="GEO">GSM847326</Accession>
#   </Sample>
# 
#   <Sample iid="GSM847327">
#     <Accession database="GEO">GSM847327</Accession>
#   </Sample>
# 
#   <Sample iid="GSM847328">
#     <Accession database="GEO">GSM847328</Accession>
#   </Sample>
# 
#   <Sample iid="GSM847329">
#     <Accession database="GEO">GSM847329</Accession>
#   </Sample>
# 
#   <Sample iid="GSM847330">
#     <Accession database="GEO">GSM847330</Accession>
#   </Sample>
# 
#   <Sample iid="GSM847331">
#     <Accession database="GEO">GSM847331</Accession>
#   </Sample>
# 
#   <Sample iid="GSM847332">
#     <Accession database="GEO">GSM847332</Accession>
#   </Sample>
# 
#   <Sample iid="GSM847333">
#     <Accession database="GEO">GSM847333</Accession>
#   </Sample>
# 
#   <Sample iid="GSM847334">
#     <Accession database="GEO">GSM847334</Accession>
#   </Sample>
# 
#   <Sample iid="GSM847335">
#     <Accession database="GEO">GSM847335</Accession>
#   </Sample>
# 
#   <Sample iid="GSM847336">
#     <Accession database="GEO">GSM847336</Accession>
#   </Sample>
# 
#   <Sample iid="GSM847337">
#     <Accession database="GEO">GSM847337</Accession>
#   </Sample>
# 
#   <Sample iid="GSM847338">
#     <Accession database="GEO">GSM847338</Accession>
#   </Sample>
# 
#   <Sample iid="GSM847339">
#     <Accession database="GEO">GSM847339</Accession>
#   </Sample>
# 
#   <Series iid="GSE34318">
#     <Status database="GEO">
#       <Submission-Date>2011-12-09</Submission-Date>
#       <Release-Date>2012-07-06</Release-Date>
#       <Last-Update-Date>2013-01-30</Last-Update-Date>
#     </Status>
#     <Title>Mapping regulatory elements using signatures of open chromatin in Arabidopsis thaliana</Title>
#     <Accession database="GEO">GSE34318</Accession>
#     <Pubmed-ID>22773751</Pubmed-ID>
#     <Summary>
# Gene expression is controlled by the complex interaction of transcription factors binding to promoters and other regulatory DNA elements.  One common characteristic of the genomic regions associated with regulatory proteins is a pronounced sensitivity to DNase I digestion.  We reported genome-wide high resolution maps of DNase I hypersensitive (DH) sites from both seedling and flower tissues of Arabidopsis from the Columbia (Col) ecotype and the corresponding ddm1 (deficient in DNA methylation 1) mutant. We identified 38,290, 41,193, 38,313, and 38,153 DH sites in leaf (Col), flower (Col), ddm1 leaf, and ddm1 flower tissues, respectively. Approximately 45% of the DH sites in all tissue types were located within 1 kb of a transcription start site (TSS), which represents a putative promoter region. Pairwise comparisons of the DH sites derived from different tissue types revealed DH sites specific to each tissue. DH sites are significantly associated with long non-coding RNAs (lncRNAs) and conserved non-coding sequences (CNSs). The binding sites of MADS-domain transcription factors AP1 and SEP3 are highly correlated with DH sites.
#     </Summary>
#     <Overall-Design>
# To map the DH sites in A. thaliana, we constructed a total of five DNase-seq libraries using leaf and flower tissues from the Columbia (Col) ecotype and a ddm1 (deficient in DNA methylation 1) mutant of Columbia.  These libraries were sequenced using the Illumina Genome Analyzer.  We obtained a total of 190 million (M) sequence reads from these libraries.  Approximately 114 M reads had a single sequence match in the A. thaliana genome
#     </Overall-Design>
#     <Type>Expression profiling by high throughput sequencing</Type>
#     <Type>Genome binding/occupancy profiling by high throughput sequencing</Type>
#     <Contributor-Ref ref="contrib2" position="1" />
#     <Contributor-Ref ref="contrib3" position="2" />
#     <Contributor-Ref ref="contrib1" position="3" />
#     <Contributor-Ref ref="contrib4" position="4" />
#     <Contact-Ref ref="contrib1" />
#     <Sample-Ref ref="GSM847326" />
#     <Sample-Ref ref="GSM847327" />
#     <Sample-Ref ref="GSM847328" />
#     <Sample-Ref ref="GSM847329" />
#     <Sample-Ref ref="GSM847330" />
#     <Sample-Ref ref="GSM847331" />
#     <Sample-Ref ref="GSM847332" />
#     <Sample-Ref ref="GSM847333" />
#     <Sample-Ref ref="GSM847334" />
#     <Sample-Ref ref="GSM847335" />
#     <Sample-Ref ref="GSM847336" />
#     <Sample-Ref ref="GSM847337" />
#     <Sample-Ref ref="GSM847338" />
#     <Sample-Ref ref="GSM847339" />
#     <Supplementary-Data type="TAR">
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE34nnn/GSE34318/suppl/GSE34318_RAW.tar
#     </Supplementary-Data>
#     <Supplementary-Data type="BED">
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE34nnn/GSE34318/suppl/GSE34318_dhsites_region.bed.gz
#     </Supplementary-Data>
#     <Supplementary-Data type="TXT">
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE34nnn/GSE34318/suppl/GSE34318_fpkm_exp.txt.gz
#     </Supplementary-Data>
#     <Supplementary-Data type="SRA Study">
# ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP009/SRP009678
#     </Supplementary-Data>
#     <Relation type="SRA" target="http://www.ncbi.nlm.nih.gov/sra?term=SRP009678" />
#     <Relation type="BioProject" target="http://www.ncbi.nlm.nih.gov/bioproject/PRJNA151473" />
#   </Series>
# 
# </MINiML>
