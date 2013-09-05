#!/usr/bin/env perl
use v5.12.0;
use warnings FATAL => "all";
use autodie;
use Data::Dumper;
use Cwd qw/getcwd/;
use File::Basename qw/basename dirname/;
use File::Path qw/make_path remove_tree/;
use File::Spec::Functions qw/rel2abs canonpath catdir catfile updir/;
use File::Copy;

use FindBin;
use lib "$FindBin::Bin/lib";
use aliased PFMLauncher => 'PFM';

use Pod::Usage;
use Getopt::Long;

my $result = GetOptions (
    "config|c=s"          => \(my $conf_file),
    "ecotype-a|ea=s"      => \(my $ecotype_a),
    "ecotype-b|eb=s"      => \(my $ecotype_b),
    "bowtie-index-a|ba=s" => \(my $bowtie_index_a),
    "bowtie-index-b|bb=s" => \(my $bowtie_index_b),
    "reads|r=s"           => \(my $reads),
    "reads-label|l=s"     => \(my $reads_label),
    "annotation|a=s"      => \(my $annotation),
    "skip-novel-detection|s" => \(my $skip_novel_detection),
);
pod2usage(-verbose => 2, -noperldoc => 1) if (
    !$result 
    || ! $conf_file      || ! -f $conf_file
    || ! $reads          || ! -f $reads
    || ! $annotation     || ! -f $annotation
    || ! $ecotype_a      || ! $ecotype_b
    || ! $bowtie_index_a || ! $bowtie_index_b
    || ! $reads_label
);  

my $output_dir_a = "$reads_label-vs-$ecotype_a";
my $output_dir_b = "$reads_label-vs-$ecotype_b";
my $full_cufflinks_dir = "$output_dir_a-cufflinks";
my $split_cufflinks_dir_a = "$output_dir_a-split-cufflinks";
my $split_cufflinks_dir_b = "$output_dir_b-split-cufflinks";
my $accepted_hits_a = catfile($output_dir_a, "accepted_hits.bam");
my $accepted_hits_b = catfile($output_dir_b, "accepted_hits.bam");
my $cuffdiff_dir = "$reads_label-cuffdiff";

my $transcriptome = "$annotation-transcriptome";

say "starting divorce with " .  PFM->processes . " processes";


#######################################################################
# run first two tophats

my $wait_in_between_tophats = ! -f "$transcriptome.fa";

PFM->launch(
    create_cmd($conf_file, 'tophat', {
            output                => $output_dir_a,
            GTF                   => $annotation,
            'transcriptome-index' => $transcriptome,
        }, [$bowtie_index_a, $reads]),
    expected => catfile($output_dir_a, "accepted_hits.bam")
);

PFM->wait_all_children if $wait_in_between_tophats;

PFM->launch(
    create_cmd($conf_file, 'tophat', {
            output                => $output_dir_b,
            GTF                   => $annotation,
            'transcriptome-index' => $transcriptome,
        }, [$bowtie_index_b, $reads]),
    expected => catfile($output_dir_b, "accepted_hits.bam")
);
PFM->wait_all_children;

#######################################################################
# split 

my $filtered_sam_a = catfile($output_dir_a, "everything.sam.filtered");
my $filtered_sam_b = catfile($output_dir_b, "everything.sam.filtered");
my $filtered_bam_a = catfile($output_dir_a, "everything.sam.filtered.bam");
my $filtered_bam_b = catfile($output_dir_b, "everything.sam.filtered.bam");

prepare_tophat_for_soms($output_dir_a);
prepare_tophat_for_soms($output_dir_b);
PFM->wait_all_children;
PFM->launch(
    join(" ", 
        "split_on_mismatches_sam.pl", 
        catfile($output_dir_a, "everything.sam"), 
        catfile($output_dir_b, "everything.sam"), 
    ),
    expected => [ $filtered_sam_a, $filtered_sam_b ],
);

PFM->wait_all_children;

PFM->launch("perl -S sam2bam.pl $filtered_sam_a", expected => $filtered_bam_a);
PFM->launch("perl -S sam2bam.pl $filtered_sam_b", expected => $filtered_bam_b);

PFM->wait_all_children;


#######################################################################
# cufflinks - run ecotype a accepted_hits.bam to create transcripts

my $annotation_plus_novel_transcripts;
if ($skip_novel_detection){
    $annotation_plus_novel_transcripts = $annotation;
}
else{
    $annotation_plus_novel_transcripts = catfile($full_cufflinks_dir, "transcripts.gtf");

    PFM->launch(
        create_cmd($conf_file, 'cufflinks', {
                'output-dir' => $full_cufflinks_dir,
                GTF => $annotation,
            }, [$accepted_hits_a],
        ),
        expected => [
            catfile($full_cufflinks_dir, "genes.fpkm_tracking"),  
            catfile($full_cufflinks_dir, "isoforms.fpkm_tracking"),  
            catfile($full_cufflinks_dir, "skipped.gtf"),  
            $annotation_plus_novel_transcripts,
        ]
    );

    PFM->wait_all_children;
}

#######################################################################
# cuffdiffs of split files against each other with respect to the 
# annotation+detected novel transcripts, from above

PFM->launch_and_wait(
    "cuffdiff -o $cuffdiff_dir --library-type fr-firststrand --library-norm-method classic-fpkm $annotation_plus_novel_transcripts $filtered_bam_a $filtered_bam_b",
    expected => [map { catfile($cuffdiff_dir, $_)} qw{bias_params.info
        cds.count_tracking cds.diff cds_exp.diff cds.fpkm_tracking
        cds.read_group_tracking gene_exp.diff genes.count_tracking
        genes.fpkm_tracking genes.read_group_tracking isoform_exp.diff
        isoforms.count_tracking isoforms.fpkm_tracking
        isoforms.read_group_tracking promoters.diff read_groups.info run.info
        splicing.diff tss_group_exp.diff tss_groups.count_tracking
        tss_groups.fpkm_tracking tss_groups.read_group_tracking var_model.info}
        ],
);

#     create_cmd($conf_file, 'cufflinks', {
#             'output-dir' => $split_cufflinks_dir_a,
#             GTF => $annotation_plus_novel_transcripts,
#         }, [$filtered_bam_a],
#     ),
#     expected => [
#         catfile($split_cufflinks_dir_a, "genes.fpkm_tracking"),  
#         catfile($split_cufflinks_dir_a, "isoforms.fpkm_tracking"),  
#         catfile($split_cufflinks_dir_a, "skipped.gtf"),  
#         catfile($split_cufflinks_dir_a, "transcripts.gtf"),
#     ],
# );

#######################################################################
# ...


# create_cmd('params.yaml', 'tophat', { x => 1 }, ['foo.fastq']
sub create_cmd{
    use YAML qw/LoadFile/;
    my ($yaml, $application, $additional_opts, $additional_bare) = @_;

    die "additional_opts, if present, should be href" if $additional_opts && ref $additional_opts ne 'HASH';
    die "additional_bare, if present, should be aref" if $additional_bare && ref $additional_bare ne 'ARRAY';

    my $config_href = LoadFile($yaml);

    if (exists $config_href->{$application}){
        my %application_hash = %{$config_href->{$application}};

        my @accum = ($application);
        my $bare_aref = delete $application_hash{__bare__};

        # value can be empty
        push @accum, map { ("--$_", $application_hash{$_} || ()) } keys %application_hash;

        if ($additional_opts){
            push @accum, map { ("--$_", $additional_opts->{$_} || ()) } keys %$additional_opts;
        }

        if ($bare_aref){
            push @accum, @$bare_aref;
        }
        if ($additional_bare){
            push @accum, @$additional_bare;
        }
        return join " ", @accum;
    }
    else{
        die "no such app $application in $yaml";
    }
}

sub prepare_tophat_for_soms{
    my $tophat_output_dir = shift;
    if (! -d $tophat_output_dir){
        die "no such directory $tophat_output_dir";
    }
    my $bam_mapped     = catfile($tophat_output_dir, "accepted_hits.bam");
    my $bam_unmapped   = catfile($tophat_output_dir, "unmapped.bam");

    my $sam_mapped     = catfile($tophat_output_dir, "accepted_hits.sam");
    my $sam_unmapped   = catfile($tophat_output_dir, "unmapped.sam");
    my $sam_everything = catfile($tophat_output_dir, "everything.sam");

    if (! -f $bam_mapped || ! -f $bam_unmapped){
        die "$bam_mapped or $bam_unmapped do not exist?";
    }
    if (-f $sam_everything){
        return;
    }

    PFM->launch("samtools view -h -o $sam_mapped $bam_mapped",     expected => $sam_mapped);
    PFM->launch("samtools view -h -o $sam_unmapped $bam_unmapped", expected => $sam_unmapped);
    PFM->wait_all_children;
    PFM->launch_and_wait("sam-combine.pl -o $sam_everything $sam_unmapped $sam_mapped", expected => $sam_everything);
    PFM->launch_and_wait("sam-sort.pl -r $sam_everything");
}

=head1 NAME

divorce_imprinting_tophat.pl - divorce imprinting with tophat instead of bowtie.

=head1 SYNOPSIS

 divorce_imprinting_tophat.pl -c params.yaml -eb kit7 -ea nip7 -bb msu7.0/Rice_7.0_KitCM -ba msu7.0/Rice_7.0_NippCM -r raw/mini.fastq -l mini -a msu7.0/all.chr_renamed.gff3

=head1 OPTIONS

=over

=item --reads <reads.fastq> | -r <reads.fastq>

=item -l <label>

=item -a <annotation>

=item --skip-novel-detection | -s

Don't run cufflinks on ecotype-a to create a new transcript including novelly detected transcripts.

=item --config <file.yaml> | -c <file.yaml>

Config file. 

=item -ea <name> | --ecotype-a <name>

=item -eb <name> | --ecotype-b <name>

=item -ba <bowtie_index_prefix> | --bowtie-index-a <bowtie_index_prefix>

=item -bb <bowtie_index_prefix> | --bowtie-index-b <bowtie_index_prefix>

=back

=cut

