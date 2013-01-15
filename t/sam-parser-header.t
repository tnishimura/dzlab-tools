#!/usr/bin/env perl
use strict;
use warnings;
use 5.010_000;
use Data::Dumper;
use autodie;
use Test::More;
use Test::Exception;
use YAML qw/Load Dump LoadFile DumpFile/;
use Sam::Parser;

{
    my $parser = Sam::Parser->new(file => \*DATA);
    is($parser->sam_version(), "1.0", "sam_version");
    is($parser->sort_order(),  "unsorted", "unsorted");
    is($parser->program_name(),  "Bowtie", "program_name");
    is($parser->program_version(),  "0.12.7", "program_version");
    is($parser->length->{CHR1}, 30432563, 'chr1');
    is($parser->length->{CHR2}, 19705359, 'chr2');
    is($parser->length->{CHR3}, 23470805, 'chr3');
    is($parser->length->{CHR4}, 18585042, 'chr4');
    is($parser->length->{CHR5}, 26992728, 'chr5');
    is($parser->length->{CHRC}, 154478, 'chrc');
    is($parser->length->{CHRM}, 366924, 'chrm');
}

# {
    # my $contents = $samfiles->{rc_length};
    # my $parser = Sam::Parser->new(file => \$contents, convert_rc => 1);
    # is($parser->{sam_version}, "1.0", "sam_version");
    # is($parser->{sort_order},  "unsorted", "unsorted");
    # is($parser->{length}{chr1}, $parser->{length}{RC_chr1}, 'chr1 forward/reverse matches');
    # is($parser->{length}{chr2}, $parser->{length}{RC_chr2}, 'chr2 forward/reverse matches');
    # is($parser->{length}{chr3}, $parser->{length}{RC_chr3}, 'chr3 forward/reverse matches');
    # is($parser->{length}{chr4}, $parser->{length}{RC_chr4}, 'chr4 forward/reverse matches');
    # is($parser->{length}{chr5}, $parser->{length}{RC_chr5}, 'chr5 forward/reverse matches');
    # is($parser->{length}{chrc}, $parser->{length}{RC_chrc}, 'chrc forward/reverse matches');
    # is($parser->{length}{chrm}, $parser->{length}{RC_chrm}, 'chrm forward/reverse matches');
    # ok(exists $parser->{rc_sequence}{RC_chr1}, "RC_chr1 exists in rc_sequence");
    # ok(exists $parser->{rc_sequence}{RC_chr2}, "RC_chr2 exists in rc_sequence");
    # ok(exists $parser->{rc_sequence}{RC_chr3}, "RC_chr3 exists in rc_sequence");
    # ok(exists $parser->{rc_sequence}{RC_chr4}, "RC_chr4 exists in rc_sequence");
    # ok(exists $parser->{rc_sequence}{RC_chr5}, "RC_chr5 exists in rc_sequence");
    # ok(exists $parser->{rc_sequence}{RC_chrc}, "RC_chrc exists in rc_sequence");
    # ok(exists $parser->{rc_sequence}{RC_chrm}, "RC_chrm exists in rc_sequence");
# }


#open my $handle_rc_length, '<', $samfiles{$samfiles->{rc_length}};
#close $handle_rc_length;

done_testing();
__DATA__
@HD	VN:1.0	SO:unsorted
@SQ	SN:chr1	LN:30432563
@SQ	SN:chr2	LN:19705359
@SQ	SN:chr3	LN:23470805
@SQ	SN:chr4	LN:18585042
@SQ	SN:chr5	LN:26992728
@SQ	SN:chrc	LN:154478
@SQ	SN:chrm	LN:366924
@SQ	SN:RC_chr1	LN:30432563
@SQ	SN:RC_chr2	LN:19705359
@SQ	SN:RC_chr3	LN:23470805
@SQ	SN:RC_chr4	LN:18585042
@SQ	SN:RC_chr5	LN:26992728
@SQ	SN:RC_chrc	LN:154478
@SQ	SN:RC_chrm	LN:366924
@PG	ID:Bowtie	VN:0.12.7	CL:"bowtie -S -f -B 1 -v 3 --best /wip/tools/genomes/AT/TAIR_reference.fas read.rc.fasta"
