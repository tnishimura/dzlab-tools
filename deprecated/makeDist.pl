#!/usr/bin/perl
# makeDist.pl
# written by Tracy Ballinger, last edited 1/15/08
# added another stop option, option 6
# 5/15/06- added a highestend to get rid of effects of overlapping genes, and added lastloci to backup if loci are missed in the range. 
# Modified by Pedro Silva (psilva@dzlab.pmb.berkeley.edu), last edited 09/16/08

use warnings;
use strict;
use List::Util qw[min max];
use Switch;

# Check for necessary arguments.
# No arguments calls usage
if (@ARGV <7) {
    die("\nUSAGE:\nmakeDist.pl <genefile> <ratiofile> <[outputfile][-]> <bin width> <distance> <stop flag> [stop distance] <3prime/5prime>\n\nTakes standard GFF files as input.\naccepts multiple chromosomes per file, and sorts them per start location if needed.\n\n");
}

##### GET PARAMETERS #####

# Read in gene transcript file (string)
# This is a standard, per spec, GFF file with 9 fields:
# <chr><source><type><lstart><lend><score><strand><phase><attributes>
# <chr>, <source>, <type>, <score>, and <phase> are discarded. This should be cleaned up in the future.
my $genefile = "$ARGV[0]";
open (GNE, "<$genefile");
my @genes=<GNE>;
close(GNE);

# Read in features data file/ratio file from nimblegen arrays (string)
# This is a standard, per spec, GFF file with 9 fields:
# <chr><source><type><lstart><lend><score><strand><phase><attributes>
# <source>, <type>, <strand>, <phase>, and <attributes> are discarded. This should be cleaned up in the future.
my $runsfile = "$ARGV[1]";
open (RNS, "<$runsfile");
@main::locs = <RNS>;
close(RNS);

# Read in name of output file (string)
# Will append .ends suffix in actual output
# the distfile ends in .ends and will list scores for each 
# individual gene transcript
my $minarg=0;
if($ARGV[2] ne "-") {
    my $outfname="$ARGV[2]";
    my $distfile = "$outfname.ends";
    open(DIS, ">$distfile");
}
else {
    $minarg=1;
    print STDERR "Printing to standard output\n";
}

# Sliding window for analysis in bp (int)
my $binwidth=$ARGV[3];

# Half-maximum amplitude of the analysis in bp (int)
my $distance = $ARGV[4];

# Analysis will go to a maximum of distance*2 (depending on flags)
# Divide 2*distance by binwidth for the total number of scores in the analysis
my $numbins = 2 * int($distance / $binwidth);

# special stop options (int)
# stop==0: don't look at adjacent genes, but just keep going the maximun distance away 
# stop==2: look within the gene, but look upstream or downstream in intergenic regions (stop at adjacent gene).  So, when aligning at the 5' end, will stop at the 3' end, or at an adjacent gene, whichever comes first, but will not include intergenic regions downstream.  It will include intergenic regions upstream of the 5' end.  
# stop==1: look at adjacent genes and stop looking only when you encounter another gene. ie. does NOT stop at the 3' end when aligning at the 5' end, but keeps going into intergenic regions until it hits another gene. 
# stop==3: only look within a gene, so stop at the end of the gene.
# stop==4: start at the start of the gene if there's an overlap, otherwise go upstream.  stop at the end of the gene. 
# stop==5: start the given distance upstream of the gene end (ignore adjacent or overlapping genes), and stop at the end of the gene.  So if aligning at the 5' end, will start x bp upstream and stop at the 3' end.
# stop==6: the same as stop flag #2, but stopping alignment by the user-specified amount before the end of the gene.
my $stop=$ARGV[5];

# $stopdistance is for flag 2 only. It is the distance before the end of a gene where alignment will stop (int)
my $stopdistance=0;

# 3prime or 5prime alignment (string)
my $prime;

if($stop==6) {
    if(@ARGV<8) {
	die "Stop flag 6 requires the stop distance away from the gene in bp as an additional argument after it\n";
    }
    $stopdistance=$ARGV[6];
    $prime="$ARGV[7]";
}
else {
    # 3prime or 5prime (string)
    $prime="$ARGV[6]";
}

# Check that alignment option if valid
if ($prime ne "5prime" && $prime ne "3prime" && $prime ne "middle"){
    die "The last argument must be either \"3prime\", \"5prime\" or \"middle\"\n";
}

# open up a log file to store analysis information
my $logfile="$ARGV[2].log";
open(LOG, ">$logfile");

print LOG "# genefile:\t", $genefile, "\n";
print LOG "# datafile:\t", $runsfile, "\n";
print LOG "# endsfile:\t", $ARGV[2], "\n";
print LOG "# binwidth:\t", $binwidth, "\n";
print LOG "# distance:\t", $distance, "\n";
print LOG "# stopflag:\t", $stop, "\n";
print LOG "# stopdistance: ", $stopdistance, "\n" unless $stop!=6;
print LOG "# alignment:\t", $prime, "\n\n";


#####PRE-PROCESSING#####

# finds rows in input arrays that have changes in chromosomes
my @genesindex=findArray(@genes);
my @locsindex=findArray(@main::locs);

print STDERR "Number of chromosomes in $genefile:\t", scalar(@genesindex), "\n";
print STDERR "Number of chromosomes in $runsfile:\t", scalar(@locsindex), "\n";
print LOG "Number of chromosomes in $genefile:\t", scalar(@genesindex), "\n";
print LOG "Number of chromosomes in $runsfile:\t", scalar(@locsindex), "\n";

# saves our genes and locs files for later recall;
my @genesbackup=@genes;
my @locsbackup=@main::locs;

# goes through all tracked change rows, clears out @tmp and
# assigns the correct range of elements from input array to @tmp
for(my $l=0;$l<@genesindex;$l++) {

    $main::lastlocindex=0;
    
    # this is in case a gene overlaps a previous gene such that the previous gene is not necessarily the closest one. 
    my $highestend=0;

    print STDERR "Found chromosome $l in $genefile at line $genesindex[$l]:\t", (split '\t', $genesbackup[$genesindex[$l]])[0], "\n";
    print STDERR "Found chromosome $l in $runsfile at line $locsindex[$l]:\t", (split '\t', $locsbackup[$locsindex[$l]])[0], "\n";
    print LOG "Found chromosome $l in $genefile at line $genesindex[$l]:\t", (split '\t', $genesbackup[$genesindex[$l]])[0], "\n";
    print LOG "Found chromosome $l in $runsfile at line $locsindex[$l]:\t", (split '\t', $locsbackup[$locsindex[$l]])[0], "\n";

    # clears and puts into @genes a sub array with all equal IDs
    if ($l+1<@genesindex) {	
	@genes=sortArray(splitArray($genesindex[$l], $genesindex[$l+1], @genesbackup));
    }
    else {
	@genes=sortArray(splitArray($genesindex[$l], scalar(@genesbackup)-1, @genesbackup));
    }
    # clears and puts into @locs a sub array with all equal IDs
    if ($l+1<@locsindex) {
	@main::locs=sortArray(splitArray($locsindex[$l], $locsindex[$l+1], @locsbackup));
    }
    else {
	@main::locs=sortArray(splitArray($locsindex[$l], scalar(@locsbackup)-1, @locsbackup));
    }


#--------Actual doing stuff program-------------------------------- 
# initialize variables for neighboring genes
    my @previous=(0,0,0,0,0);
    my @next =(0,0,0,0,0);
    my ($prime3, $prime5);
    my $localmin=0;
    my $localmax=0;
    for(my $j=0;$j<@genes; $j++){
	# 2008/08/14 note: this script now expects the following gene file form:
	# <chr><source><type><lstart><lend><score><strand><phase><attributes>
	# <chr>, <source>, <type>, <score>, and <phase> are discarded.
	my @tempparameters = split(/\t|\r/, $genes[$j]);
	our ($start, $end, $strand, $id) = (@tempparameters[3,4,6,8]);

	if ( ($prime eq "5prime" && $strand eq "-") || ($prime eq "3prime" && $strand eq "+") ) {
	    # switch to 3prime analysis
	    $prime5 = $end;
	    $prime3 = $start;
	}
	else {
	    # switch to 5prime analysis
	    $prime5 = $start;
	    $prime3 = $end;
	}
	
	if ($prime eq "middle") {
	    # switch to middle analysis
	    $prime5 = int(($end-$start+1)/2) + $start;
	}
	
	my $rangestart = $prime5 - $distance;
	my $rangeend = $prime5 + $distance;
	my ($pend, $nstart);
	
	# use the stop flag for stopping at adjacent genes
	switch ($stop) {
	    case /[126]/ {
		# initialize variables for neighboring genes
		@previous=split(/\t/,$genes[$j-1]) unless ($j==0);
		if ($previous[4]>$highestend){
		    $pend=$previous[4];
		    $highestend=$pend;
		}
		else{$pend=$highestend;}
		if ($j==$#genes){
		    $nstart=0;
		}
		else {
		    @next=split(/\t/,$genes[$j+1]);
		    $nstart=$next[3];}
		if ($pend>$rangestart){ # pend is always >=0
		    $rangestart=$pend;}
		if ($nstart<$rangeend && $nstart !=0){
		    $rangeend=$nstart;}
		# for stop==2, still stop at adjacent genes, but also stop at the end of the gene
		if ($stop==2){
		    if (($prime eq "5prime" && $strand eq "+") || ($prime eq "3prime" && $strand eq "-")){ #assigning the end of the range
			if ($end < $rangeend){
			    $rangeend=$end}
		    }
		    elsif(($prime eq "3prime" && $strand eq "+") || ($prime eq "5prime" && $strand eq "-")) {
			if ($start > $rangestart) {
			    $rangestart=$start;
			}
		    }
		    elsif ($prime eq "middle"){
			die "Cannot pick middle with stop option 2\n";
		    }
		}
		# for stop==6, still stop at adjacent genes, but also stop at the end of the gene minus the user-specified stop distance in bp
		if ($stop==6) {
		    if (($prime eq "5prime" && $strand eq "+") || ($prime eq "3prime" && $strand eq "-")){ #assigning the end of the range
			if ($end-$stopdistance < $rangeend){
			    $rangeend=$end-$stopdistance;}
		    }
		    elsif(($prime eq "3prime" && $strand eq "+") || ($prime eq "5prime" && $strand eq "-")) {
			if ($start+$stopdistance > $rangestart) {
			    $rangestart=$start+$stopdistance;
			}
		    }
		    elsif ($prime eq "middle"){
			die "Cannot pick middle with stop option 6\n";
		    }
		}
	    }
	    
	    case 3 {
		if ($start>$rangestart){$rangestart=$start;}
		if ($end<$rangeend){$rangeend=$end;}
	    }
	    
	    case 4 {
		@previous=split(/\t/,$genes[$j-1]) unless ($j==0);
		@next=split(/\t/,$genes[$j+1]) unless ($j==$#genes);
		if ($previous[4]>$highestend){
		    $pend=$previous[4];
		    $highestend=$pend;
		}
		else{$pend=$highestend;}
		$nstart=$next[0];
		if (($prime eq "5prime" && $strand eq "+") || ($prime eq "3prime" && $strand eq "-")){
		    if ($start<$pend){
			$rangestart=$start;}
		    elsif ($pend > $rangestart){ # start > pen
			$rangestart=$pend;
		    }
		    if ($end < $rangeend){
			$rangeend=$end;
		    }
		}
		elsif(($prime eq "3prime" && $strand eq "+") || ($prime eq "5prime" && $strand eq "-")){
		    if ($end > $nstart){
			$rangeend=$end;}
		    elsif ($nstart < $rangeend && $nstart != 0){
			$rangeend=$nstart;
		    }
		    if ($start > $rangestart){
			$rangestart=$start;}
		}
		elsif($prime eq "middle"){
		    die "Cannot pick middle with stop option 4\n";
		}
	    }
	    
	    case 5 {
		# want to go ignore adjacent genes upstream, and stop at the end of the gene
		if (($prime eq "5prime" && $strand eq "+") || ($prime eq "3prime" && $strand eq "-")){
		    if ($end < $rangeend){$rangeend=$end;}
		}
		elsif(($prime eq "3prime" && $strand eq "+") || ($prime eq "5prime" && $strand eq "-")){
		    if ($start > $rangestart){$rangestart=$start;}
		}
	    }
	    
	    else {
		die "Unknown stop flag.\n";
	    }
	}
	
	our ($binstart, $binend, @binstarts);
	for(my $i=0; $i<$numbins; $i++){
	    #numbins is always an even number
	    $binstart = $prime5 - (($numbins / 2)-$i) * $binwidth; 
	    $binend = $binstart + $binwidth -1;
	    if ($rangestart < $binend && $rangeend > $binstart){
		$binstarts[$i]=$binstart; 
	    }
	    else{
		$binstarts[$i]="na";
	    }
	}
	#end of stop stuff
	
	my @scorelist = &getScores(@binstarts);

        # now print stuff out 
	my $tmpwrd = (split ';', "$id")[0]; #print only until first field delimiter ' ; '
	$tmpwrd =~ s/"//g; #delete quotes ' " '
	if ($strand eq "-"){
	    my @tmpscores= reverse @scorelist;
	    @scorelist=@tmpscores;
	}
	foreach my $scr (@scorelist){
	    if ($scr eq "na"){
		$tmpwrd.="\tna";
	    }
	    else{
		$tmpwrd.= sprintf "\t%.5g", ($scr);
		$localmin=$scr if $scr<$localmin;
		$localmax=$scr if $scr>$localmax;
	    }
	}
	$tmpwrd.="\n";
	if($minarg==0) {
	    print DIS $tmpwrd; 
	}
	else {print $tmpwrd;}
    }
    print STDERR "Local minimum for ", (split '\t', $genesbackup[$genesindex[$l]])[0], ": ", $localmin, "\n";
    print STDERR "Local maximum for ", (split '\t', $genesbackup[$genesindex[$l]])[0], ": ", $localmax, "\n";
    print LOG "Local minimum for ", (split '\t', $genesbackup[$genesindex[$l]])[0], ": ", $localmin, "\n";
    print LOG "Local maximum for ", (split '\t', $genesbackup[$genesindex[$l]])[0], ": ", $localmax, "\n";
}

close(LOG);
close(DIS);

exit (0);

#-----------Find indices for each ID change---------------#
# finds each index that signifies a change in type of record
# returns index array
sub findArray {
    my @array=@_;

    # @index contains a list of locations where array should be split
    # $previousid and $currentid are scalars containing the previous and current IDs (chr1, 2, etc)
    # $chrcount is just a counter for the number of different types of IDs
    my (@index, $previousid, $currentid, @tmp);
    my $chrcount=0;
    
    # goes through full gene file
    for (my $k=0;$k<@array;$k++) {

	# gets current starting coordinate
	$currentid=(split '\t', $array[$k])[0]; #chr1, chr2, etc
	
	# if we're at beginning of file if doesn't make sense to look for changes already
	# gets previous starting coordinate
	if ($k!=0) {$previousid=(split '\t', $array[$k-1])[0];}	
	else {$previousid=$currentid;}
	
	# keeps track of number of different types of records
	# also stores each record type change in @index
	# ignores pound (#) characters if they're the first printing character in the record
	if ( ($currentid ne $previousid && $currentid !~ m-^\s*#-) || $k==0 )  {
	    if($currentid !~ m-^\s*#-) {
		$index[$chrcount]=$k;
		$chrcount++;
	    }
	}
    }
    return @index;
}


#----------Split into multiple arrays----------------#
# takes an input index array
# returns split array
sub splitArray {
    my ($start, $end)=($_[0], $_[1]);
    my @array=@_[2..@_];
    return @array[$start..$end-1];
}


#-----------------sortArray-------------------------------
# sorts input array at fourth field delimited by tabs (starting coordinates)
# returns array sorted numerically on the 4th field
sub sortArray {
    return sort {
	(split '\t', $a)[3] <=> (split '\t', $b)[3]
    } @_;
}


#-----------------getScores--------------------------------
sub getScores {
    # rangestart <= prime <= rangeend
    my @binstarts=@_;
    my @binsums = my @binscores = my @binfinal=();    
    my $firsti = my $firstloci = -1; 
    
    for(my $j=0; $j<@binstarts; $j++){
        my $xstart = $binstarts[$j];
        $binsums[$j]=0;
        $binscores[$j]=0;
        # if it's na, don't look
        if ($xstart eq "na"){
            $binsums[$j]="na";
            $binscores[$j]="na"; 
            next;
        }
        my $xend = $xstart + $binwidth -1;
        my $i=$main::lastlocindex;
	
	while ($i < @main::locs){
	    my @tempparameters = split(/\t|\r|\n/, $main::locs[$i]);
            my ($lstart, $lend, $score) = (@tempparameters[3,4,5]);
	    my ($lo, $hi, $overlap, $frac);
	    
            #we've started past the first bin
            if ($firstloci==-1 && $lstart >$xstart && $i>0){
                #need to back up the probes
                $i--;
                next;
            }
            else{
                $firstloci=1;}
            
            # see if the loc overlaps the range at all
            if ($lend >= $xstart && $lstart <= $xend){ 
		
                # keep track of the index of the first loci in the whole gene
                if ($firsti == -1){
                    $firsti=$i;}
                if ($lstart < $xstart){
                    $lo = $xstart;}
                else {
                    $lo = $lstart;}
                if ($lend < $xend){
                    $hi = $lend;}
                else {
                    $hi = $xend;}
                $overlap = $hi -$lo +1; 
                $frac = $overlap / ($lend - $lstart + 1); 
                $binsums[$j] +=$frac; 
                $binscores[$j] += $frac * $score;
            }
            elsif ($lstart > $xend){
                last;
            }
            # get here if $lend < $xend ie you are still looking for 
            # runs in the range of the 5prime end
            $i++;
        } # end of while loop
        
	if ($firsti >= 0){
            $main::lastlocindex=$firsti;
        } 
    } # end of foreach binstarts
    
    for(my $i=0; $i<@binsums; $i++){
        if ($binsums[$i] eq "0" || $binsums[$i] eq "na"){
            $binfinal[$i]="na";
        }
        else {
            $binfinal[$i]=$binscores[$i] / $binsums[$i];  
        }
    }
    return @binfinal;
}
