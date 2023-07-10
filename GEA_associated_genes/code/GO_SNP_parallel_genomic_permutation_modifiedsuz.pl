#!/usr/bin/perl -w
use strict;
use Getopt::Std;

###

my $SNPFile = 'snp_set_100M.txt';

###

#Get command line arguments
our($opt_i,$opt_r, $opt_b, $opt_e);
getopt('irbe');

unless (defined($opt_r) && defined($opt_i) && defined($opt_e)){
    print <<EOF;

SNP-level GO enrichment analysis

It requires the following arguments:
-i [input file prefix]
-r [number of replicates to simulate]
-e [path prefix of output gene list for each environmental variable]
-b [OPTIONAL - batch number, optional to allow parallel processing]

There must be an input file for each chromosome arm with the above prefix followed by "_Chr2L.txt", etc.
Thos files should have 1 header row, outlier site numbers in the first column (col 0), and outlier statistic values in the second column (col 1).
Output files will have the above prefix appended by "_GO[BatchNum].txt".

EOF
exit 1;
}

#Get command line arguments
my $OutlierFilePrefix = "$opt_i";
my $PermutationReps = "$opt_r";
my $env_var = "$opt_e";

my $BatchNum = "";
if (defined($opt_b)){
  $BatchNum = $opt_b;
}

#GO enrichment analysis using all genes overlapping each outlier region center

#my $InputHeaderRows = 1;  #how many header rows are present in the input file before the data

my $OutlierDistanceThreshold = 20000;

my $OutliersHigh = -1;  #Set to 1 if higher values of your statistic are more extreme, or -1 if lower values are more extreme

my $ExonFilePrefix = '/home/siyuan/jobs/suzukii_WGS/EAA/downstream/GSEA/genomic_annot/anno_'; #assumes full name of annotation file is something like annot_chr2L.txt

my $GeneGOCatFile = '/home/siyuan/jobs/suzukii_WGS/EAA/downstream/GSEA/GO_annot/dmel/GO_gene_cats_parents20210708_sorted.txt';

my $GOCatDescFile = '/home/siyuan/jobs/suzukii_WGS/EAA/downstream/GSEA/GO_annot/dmel/GO_desc12_condensed_parents.txt';

#my $HighestFbgn = 265624;

my $OutputFile = $env_var . '_gene_list.txt'; # list of genes associated with all trimmed outliers
my $Output_outlier_gene = $env_var . '_SNP_gene_table.txt'; # list of trimmed outlier SNPs and associated genes
#my $PermutationReps = 100000;

my @chrs = ('Chr2L','Chr2R','Chr3L','Chr3R','ChrX');  #keep in alphatbetically sorted order
#my @StartSites = (500000,5200000,600000,6900000,2300000);  #high recomb boundaries
#my @StopSites = (17500000,20800000,17700000,26600000,21400000);
my @StartSites = (1,1,1,1,1);  #no recomb cutoffs
my @StopSites = (23011543,21146707,24543556,27905052,22422826);  #0-based
die if (@chrs != @StartSites);
die if (@chrs != @StopSites);

my $c = 0;
my $f = 0;
my $i = 0;
my $s = 0;
my $file = '';
my $center = 0;
my @line = ();
my @ChrList = ();
#my @PosList = ();
my @CDSAoA = ();
my @AllOverlapGenes = ();
my @AllOverlapFbgns = ();
my @AllWinStarts = ();
my @AllWinStops = ();
my @AllOverlapFbgnAoA = ();
my @GeneOrientationList = ();
my @GenePositionList = ();
my $j = 0;
my $k = 0;
my $g = 0;
my $h = 0;
#my $NearestExon = -1;
my $NearestDistance = 1000000000;
my $distance = 0;
#my @GenesWithin = ();
my @UpstreamGenes = ();
my @DownstreamGenes = ();
my @OverlapGenes = ();
my $AlreadyListed = 0;
my $orientation = "";
my $GenePosition = 0;
my $UpstreamLimit = 0;
my $DownstreamLimit = 0;
my $direction = 0;
my $GeneList = "";
my @OverlapFbgns = ();
my $FbgnList = "";
my $ClosestGene = "";
my $start = 0;


#For each chromosome arm, read in all coding sequence blocks
my $LastPos = 0;
my @CdsFbgns = ();
my @CdsGenes = ();
my @CdsStarts = ();
my @CdsStops = ();
my @CdsFbgnAoA = ();
my @CdsGeneAoA = ();
my @CdsStartAoA = ();
my @CdsStopAoA = ();
my @OutlierSites = ();
my @OutlierValues = ();
my @OutlierSiteAoA = ();
for ($c = 0; $c < @chrs; $c++){
  @CdsFbgns = ();
  @CdsGenes = ();
  @CdsStarts = ();
  @CdsStops = ();
  $file = $ExonFilePrefix . $chrs[$c] . '.txt';
  open E, "<$file" or die "can not open exon file $file\n";
  while (<E>){
    chomp;
    last if m/^$/;
    @line = split;
#    $line[2] =~ s/c/C/;
    push @CdsFbgns, $line[0];
    push @CdsGenes, $line[1];
    push @CdsStarts, $line[4];
    push @CdsStops, $line[5];
  }
  close E;
#Trim lists to only the regions analyzed on this arm, then add this arm's arrays to AoAs
  for ($i = 0; $i < @CdsStarts; $i++){
    if ($CdsStarts[$i] < $StartSites[$c]){
      $i = $i - 1;
      splice @CdsFbgns, 0, $i;
      splice @CdsGenes, 0, $i;
      splice @CdsStarts, 0, $i;
      splice @CdsStops, 0, $i;
      last;
    }
  }
  for ($i = 0; $i < @CdsStarts; $i++){
    if ($CdsStops[$i] > $StopSites[$c]){
      $i++;
      $j = @CdsStarts - $i;
      splice @CdsFbgns, $i, $j;
      splice @CdsGenes, $i, $j;
      splice @CdsStarts, $i, $j;
      splice @CdsStops, $i, $j;
      last;
    }
  }
  push @CdsFbgnAoA, [ @CdsFbgns ];
  push @CdsGeneAoA, [ @CdsGenes ];
  push @CdsStartAoA, [ @CdsStarts ];
  push @CdsStopAoA, [ @CdsStops ];

#Get outlier site positions and statistic values
  @OutlierSites = ();
  @OutlierValues = ();
  $file = $OutlierFilePrefix . $chrs[$c] . '.txt';
  open I, "<$file" or die "can not open exon file $file\n";
  scalar (<I>);
  while (<I>){
    chomp;
    last if m/^$/;
    @line = split;
#add "next if" statements here if you need to exclude low rec sites that exist in the output file
    push @OutlierSites, $line[0];
    push @OutlierValues, $line[1];
  }
  close I;

  # # for debugging only
  # my $num_OutlierSites = @OutlierSites;
  # print "# OutlierSites of $chrs[$c]: $num_OutlierSites\n";

#Pare down closely linked outliers - if within a threshold distance, only keep the one with a more extreme value
  for ($i = 0; $i < @OutlierSites; $i++){
    $OutlierValues[$i] = $OutlierValues[$i] * $OutliersHigh;
  } 
  # $LastPos = $OutlierSites[0];
  for ($i = 1; $i < @OutlierSites; $i++){
    $LastPos = $OutlierSites[$i-1];
    if (($OutlierSites[$i] - $LastPos) < $OutlierDistanceThreshold){
      if ($OutlierValues[$i] > $OutlierValues[$i-1]){
	$i--;
	splice @OutlierSites, $i, 1;
	splice @OutlierValues, $i, 1;
      }
      else{
	splice @OutlierSites, $i, 1;
	splice @OutlierValues, $i, 1;
	$i--;
      }
    }
    # else{
    #   $LastPos = $OutlierSites[$i];
    # }
  }
  push @OutlierSiteAoA, [ @OutlierSites ];
}
# # for debugging only
# my $num_OutlierSiteAoA = @OutlierSiteAoA;
# print "# OutlierSiteAoA: $num_OutlierSiteAoA\n";

my $TotalOutliers = 0;
for ($i = 0; $i < @OutlierSiteAoA; $i++){
  $TotalOutliers += @{$OutlierSiteAoA[$i]};
}
print "Thinned outliers:  $TotalOutliers\n";

#Read GO category list into pre-AoA
my @GOAoA = ();
my @blank = ();
my $Fbgn = '';
open G, "<$GeneGOCatFile" or die "Can't open $GeneGOCatFile\n";
#$i = -1;
while (<G>){
  chomp;
  last if m/^$/;
  @line = split;
  $Fbgn = $line[0];
  $Fbgn =~ s/FBgn//;
  $Fbgn += 0;
  while ($Fbgn >= @GOAoA){    ###
    push @GOAoA, [ @blank ];  ###
  }                           ###  
#  while ($Fbgn > ($i + 1.5)){  #Current code assumes Fbgns are presented in numerical order
#    push @GOAoA, [ @blank ];
#    $i++;
#  }
  shift @line;
  shift @line;
  @{$GOAoA[$Fbgn]} = @line;   ###
#  push @GOAoA, [ @line ];
#  $i = $Fbgn;
}
close G;

# # for debugging only
# for ($i = 1; $i < @GOAoA; $i++){
#   print "# GOAoA list $i: @{$GOAoA[$i]}\n";
# }

#if (@GOAoA != ($HighestFbgn + 1)){
#  $i = @GOAoA;
#  die "Highest Fbgn is $HighestFbgn, but GOAoA has $i rows\n";
#}

#Make arrays with all GO terms (first with repeats, then two more arrays with unique terms and their corresponding counts)
my @AllGORepeats = ();
my @AllGOUnique = ();
my @AllGOCounts = ();
for ($i = 0; $i < @GOAoA; $i++){
  for ($j = 0; $j <  @{$GOAoA[$i]}; $j++){
    push @AllGORepeats, $GOAoA[$i][$j];
  }
}
@AllGORepeats = sort{ lc($a) cmp lc($b) } @AllGORepeats;
$i = @AllGORepeats;
print "Set of all genes has a total of $i GO listings,\n";

push @AllGOUnique, $AllGORepeats[0];
push @AllGOCounts, 1;
for ($i = 1; $i < @AllGORepeats; $i++){
  if ($AllGORepeats[$i] eq $AllGORepeats[$i-1]){
    $AllGOCounts[-1]++;
  }
  else{
    push @AllGOUnique, $AllGORepeats[$i];
    push @AllGOCounts, 1;
  }
}
$i = @AllGOUnique;
print "Retrieved $i distinct GO categories corresponding to Fbgns in analyzed regions.\n";

#Identify outlier Fbgns
my $overlap = 0;
my $LastCds = 0;
my $ClosestDist = 100000000;
my $ClosestFbgn = '';
my @OutlierFbgns = ();
my @OutlierFbgnAoA = ();
open O, ">$Output_outlier_gene"; # write contig, SNP position, associated genes into a table
for ($c = 0; $c < @chrs; $c++){
  $LastCds = 0;
  for ($s = 0; $s < @{$OutlierSiteAoA[$c]}; $s++){
    @OutlierFbgns = ();
    for ($i = $LastCds; $i < @{$CdsStartAoA[$c]}; $i++){
      if (($CdsStartAoA[$c][$i] > $OutlierSiteAoA[$c][$s]) || ($i == (@{$CdsStartAoA[$c]} - 1))){
	$LastCds = $i;
#go back 200 CDS and see if any overlap this site	
	$overlap = 0;
	for ($j = $i - 1; $j > ($i - 201); $j--){
	  last if ($j < 0);
	  if (($CdsStartAoA[$c][$j] < $OutlierSiteAoA[$c][$s]) && ($CdsStopAoA[$c][$j] > $OutlierSiteAoA[$c][$s])){
	    $overlap = 1;
	    if ((@OutlierFbgns == 0) || ($OutlierFbgns[-1] ne $CdsFbgnAoA[$c][$j])){
	      push @OutlierFbgns, $CdsFbgnAoA[$c][$j];
	    }
	  }
	}
#if no overlap, get the Fbgns for the nearest exons to the left (search previous 200) and right (number $i) 
	if ($overlap == 0){
	  push @OutlierFbgns, $CdsFbgnAoA[$c][$i];
	  $ClosestDist = 1000000000;
	  $ClosestFbgn = '';
	  if ($i < (@{$CdsStartAoA[$c]} - 1)){
	    for ($j = $i - 1; $j > ($i - 201); $j--){
	      last if ($j < 0);
	      if (($OutlierSiteAoA[$c][$s] - $CdsStopAoA[$c][$j]) < $ClosestDist){
		$ClosestDist = $OutlierSiteAoA[$c][$s] - $CdsStopAoA[$c][$j];
		$ClosestFbgn = $CdsFbgnAoA[$c][$j];
	      }
	    }
	    if ((@OutlierFbgns == 0) || ($OutlierFbgns[-1] ne $ClosestFbgn)){
	      if ($ClosestFbgn){
		push @OutlierFbgns, $ClosestFbgn;
	      }
	    }
	  }
	}
	last;
      }
    }
#remove duplicates from OutlierFbgns, then add to AoA
    @OutlierFbgns = sort{ lc($a) cmp lc($b) } @OutlierFbgns;
    for ($i = 1; $i < @OutlierFbgns; $i++){
      if ($OutlierFbgns[$i] eq $OutlierFbgns[$i-1]){
	splice @OutlierFbgns, $i, 1;
	$i--;
      }
    }
    push @OutlierFbgnAoA, [ @OutlierFbgns ];
    print O $chrs[$c] . "\t" . $OutlierSiteAoA[$c][$s] . "\t" . "@OutlierFbgns" . "\n";
  }
}
close O;
	 
#Eliminate duplicate genes in OutlierFbgnAoA (remove from the outlier region containing more genes)
my $match = 0;
my $l = 0;
for ($i = 0; $i < @OutlierFbgnAoA; $i++){
  for ($j = 0; $j < @{$OutlierFbgnAoA[$i]}; $j++){
    $match = 0;
    for ($k = $i + 1; $k < @OutlierFbgnAoA; $k++){
      for ($l = 0; $l < @{$OutlierFbgnAoA[$k]}; $l++){
	if ($OutlierFbgnAoA[$i][$j] eq $OutlierFbgnAoA[$k][$l]){
	  if (@{$OutlierFbgnAoA[$i]} < @{$OutlierFbgnAoA[$k]}){
	    @line = @{$OutlierFbgnAoA[$k]};
	    splice @line, $l, 1;
	    @{$OutlierFbgnAoA[$k]} = @line;
	    last;
	  }
	  else{
	    @line = @{$OutlierFbgnAoA[$i]};
	    splice @line, $j, 1;
	    @{$OutlierFbgnAoA[$i]} = @line;
	    $match = 1;
	    last;
	  }
	}
      }
      if ($match == 1){
	$j--;
	last;
      }
    }
  }
}
$k = 0;
for ($i = 0; $i < @OutlierFbgnAoA; $i++){
  $k += @{$OutlierFbgnAoA[$i]};
}
print "Obtained $k unique genes corresponding to outlier SNPs\n";

open O, ">$OutputFile";
for ($i = 0; $i < @OutlierFbgnAoA; $i++){
  for ($j = 0; $j < @{$OutlierFbgnAoA[$i]}; $j++){
    print O $OutlierFbgnAoA[$i][$j] . "\n";
  }
}
close O;