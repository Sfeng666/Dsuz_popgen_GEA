#!/usr/bin/perl -w
use strict;
use Getopt::Std;

###

my $SNPFile = 'snp_set_100M.txt';

###

#Get command line arguments
our($opt_i,$opt_r, $opt_b);
getopt('irb');

unless (defined($opt_r) && defined($opt_i)){
    print <<EOF;

SNP-level GO enrichment analysis

It requires the following arguments:
-i [input file prefix]
-r [number of replicates to simulate]
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

my $BatchNum = "";
if (defined($opt_b)){
  $BatchNum = $opt_b;
}

#GO enrichment analysis using all genes overlapping each outlier region center

#my $InputHeaderRows = 1;  #how many header rows are present in the input file before the data

my $OutlierDistanceThreshold = 20000;

my $OutliersHigh = -1;  #Set to 1 if higher values of your statistic are more extreme, or -1 if lower values are more extreme

my $ExonFilePrefix = 'anno_'; #assumes full name of annotation file is something like annot_chr2L.txt

my $GeneGOCatFile = 'GO_gene_cats_parents20210708_sorted.txt';

my $GOCatDescFile = 'GO_desc12_condensed_parents.txt';

#my $HighestFbgn = 265624;

my $OutputFile = $OutlierFilePrefix;
$OutputFile = $OutputFile . 'GO' . $BatchNum . '.txt';
#my $OutputFile = 'DGRP_clustercounts_GO_enrich_7_FSTstrict.txt';

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
  }
}
	 
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

# # for debugging only
# my $num_OutlierFbgnAoA = @OutlierFbgnAoA;
# print "# OutlierFbgnAoA: $num_OutlierFbgnAoA\n";

#For each region, obtain the list of GO categories associated with overlapping genes.  
#Eliminate duplicate GO categories in same outlier region (avoid false positive GO enrichment results due to clusters of functionally related genes)
my @RegionGOCats = ();
my @OutlierGOList = ();
my @OutlierGOCounts = ();
for ($i = 0; $i < @OutlierFbgnAoA; $i++){
  @RegionGOCats = ();

  # # for debugging only
  # print "# OutlierFbgnAoA of outlier $i: @{$OutlierFbgnAoA[$i]}\n";

  for ($j = 0; $j < @{$OutlierFbgnAoA[$i]}; $j++){
    $Fbgn = $OutlierFbgnAoA[$i][$j];
    $Fbgn =~ s/FBgn//;
    $Fbgn += 0;
    for ($k = 0; $k < @{$GOAoA[$Fbgn]}; $k++){
      push @RegionGOCats, $GOAoA[$Fbgn][$k];
    }
  }
  @RegionGOCats  = sort{ lc($a) cmp lc($b) } @RegionGOCats;

  # # for debugging only
  # my $num_RegionGOCats = @RegionGOCats;
  # print "# RegionGOCats of outlier $i: $num_RegionGOCats\n";

  for ($j = 1; $j < @RegionGOCats; $j++){
    if ($RegionGOCats[$j] eq $RegionGOCats[$j-1]){
      splice @RegionGOCats, $j, 1;
      $j--;
    }
  }
  push @OutlierGOList, @RegionGOCats;
#  push @OutlierGOMatrix, [ @RegionGOCats ];
}
$i = @OutlierGOList;
print "These are associated with a total of $i GO listings\n";
#Count how many times each GO category is hit by outliers
for ($i = 0; $i < @AllGOUnique; $i++){
  push @OutlierGOCounts, 0;
}
@OutlierGOList = sort{ lc($a) cmp lc($b) } @OutlierGOList;    
$j = 0;
for ($i = 0; $i < @OutlierGOList; $i++){
  while ($OutlierGOList[$i] ne $AllGOUnique[$j]){
    $j++;
  }
  $OutlierGOCounts[$j]++;
}
print "Empirical GO category counts established.\n";

# # for debugging only
# my $num_OutlierGOCounts = @OutlierGOCounts;
# print "# OutlierGOCounts: $num_OutlierGOCounts\n";


#PERMUTATION - Sample random outlier regions for p values
my $r = 0;
my $random = 0;
my $ArmLength = 0;
my $RandomSites = 0;
my $NewSite = 0;
my $NewContig = 0;
my @RandomSites = ();
my @RandomSiteAoA = ();
my @SampledFbgns = ();
my @SampledFbgnAoA = ();
my @SampledGOList = ();
my @SampledGOCounts = ();
my @PValues = ();
my @SummedArmLengths = ();
my @InputSites = ();
my @InLine = ();
my @InputContigs = ();

for ($i = 0; $i < @AllGOUnique; $i++){
  push @PValues, 0;
}
for ($c = 0; $c < @StartSites; $c++){
  $ArmLength += $StopSites[$c] - $StartSites[$c] + 1;
  push @SummedArmLengths, $ArmLength;
}
print "Beginning genome-wide permutations of $TotalOutliers outlier SNPs";

open S, "<$SNPFile" or die "can not open exon file $SNPFile\n";

for ($r = 0; $r < $PermutationReps; $r++){
  $RandomSites = 0;
  @SampledFbgnAoA = ();
  @SampledGOList = ();
  @SampledGOCounts = ();
  @RandomSiteAoA = ();
  for ($c = 0; $c < @OutlierSiteAoA; $c++){
    push @RandomSiteAoA, [ @blank ];
  }
  # # for debugging only
  # my $num_RandomSiteAoA = @RandomSiteAoA;
  # print "# RandomSiteAoA: $num_RandomSiteAoA\n";

#Permute outlier SNP locations. Randomly select the same number of random outlier sites as exist in the empirical data.
  while ($RandomSites < $TotalOutliers){
    # $random = int(rand($SummedArmLengths[-1]));

#Decide which chromosome this random SNP is on.     
  #   for ($c = 0; $c < @SummedArmLengths; $c++){  
  #     if ($random < $SummedArmLengths[$c]){
	# if ($c > 0){
	#   $random = $random - $SummedArmLengths[$c-1];
	# }
	# $random = $random + $StartSites[$c];
	# die if ($random < $StartSites[$c]);
	# die if ($random > $StopSites[$c]);
#Get batches of 10000 random SNPs at a time
    if (@InputSites == 0){
      for ($i = 0; $i < 10000; $i++){
	$_ = <S>;
	chomp;
	@InLine = split;
	if (@InLine == 0){
	  die "failed to retrieve data from random SNP file\n";
	}
	push @InputContigs, $InLine[0];
	push @InputSites, $InLine[1];
      }
    }

#Draw the next contig and site number
    $NewContig = shift @InputContigs;
    $NewSite =  shift @InputSites;
    for ($c = 0; $c < @chrs; $c++){
      last if ($chrs[$c] eq $NewContig);
      if ($c == (@chrs - 1)){
	die "did not find $NewContig in the @chrs array\n";
      }
    }
    die if ($NewSite < $StartSites[$c]);
    die if ($NewSite > $StopSites[$c]);
#Make sure it's not too close to another outlier SNP.
    $match = 0;
    for ($j = 0; $j < @{$RandomSiteAoA[$c]}; $j++){
        if ((abs($RandomSiteAoA[$c][$j] - $NewSite)) < $OutlierDistanceThreshold){
  $match = 1;
  last;
      }
    }
    if ($match == 0){
      push @{$RandomSiteAoA[$c]}, $NewSite;
      $RandomSites++;
    }
  }
#Then sort the random SNP positions on each chromosome
  for ($c = 0; $c < @RandomSiteAoA; $c++){
    @{$RandomSiteAoA[$c]} = sort { $a <=> $b } @{$RandomSiteAoA[$c]};
  }
###
  print "Distribution of permuted SNPs among chromosomes: ";
  for ($c = 0; $c < @RandomSiteAoA; $c++){
    $RandomSites = @{$RandomSiteAoA[$c]};
    print "$RandomSites ";
  }
  print "\n";
###
  
#Build AoA of sampled FBGNs
  for ($c = 0; $c < @chrs; $c++){
    $LastCds = 0;
    for ($s = 0; $s < @{$RandomSiteAoA[$c]}; $s++){
      @SampledFbgns = ();
      for ($i = $LastCds; $i < @{$CdsStartAoA[$c]}; $i++){
	if (($CdsStartAoA[$c][$i] > $RandomSiteAoA[$c][$s]) || ($i == (@{$CdsStartAoA[$c]} - 1))){
	  $LastCds = $i;
#go back 200 CDS and see if any overlap this site	
	  $overlap = 0;
	  for ($j = $i - 1; $j > ($i - 201); $j--){
	    last if ($j < 0);
	    if (($CdsStartAoA[$c][$j] < $RandomSiteAoA[$c][$s]) && ($CdsStopAoA[$c][$j] > $RandomSiteAoA[$c][$s])){
	      $overlap = 1;
	      if ((@SampledFbgns == 0) || ($SampledFbgns[-1] ne $CdsFbgnAoA[$c][$j])){
		push @SampledFbgns, $CdsFbgnAoA[$c][$j];
	      }
	    }
	  }
#if no overlap, get the Fbgns for the nearest exons to the left (search previous 200) and right (number $i) 
	  if ($overlap == 0){
	    push @SampledFbgns, $CdsFbgnAoA[$c][$i];
	    $ClosestDist = 1000000000;
	    $ClosestFbgn = '';
	    if ($i < (@{$CdsStartAoA[$c]} - 1)){
	      for ($j = $i - 1; $j > ($i - 201); $j--){
		last if ($j < 0);
		if (($RandomSiteAoA[$c][$s] - $CdsStopAoA[$c][$j]) < $ClosestDist){
		  $ClosestDist = $RandomSiteAoA[$c][$s] - $CdsStopAoA[$c][$j];
		  $ClosestFbgn = $CdsFbgnAoA[$c][$j];
		}
	      }
	      if ((@SampledFbgns == 0) || ($SampledFbgns[-1] ne $ClosestFbgn)){
		if ($ClosestFbgn){
		  push @SampledFbgns, $ClosestFbgn;
		}
	      }
	    }
	  }
	  last;
	}
      }
#remove duplicates from SampleFbgns, then add to AoA
      @SampledFbgns = sort{ lc($a) cmp lc($b) } @SampledFbgns;
      for ($i = 1; $i < @SampledFbgns; $i++){
	if ($SampledFbgns[$i] eq $SampledFbgns[$i-1]){
	  splice @SampledFbgns, $i, 1;
	  $i--;
	}
      }
      push @SampledFbgnAoA, [ @SampledFbgns ];
    }
  }
#  $k = 0;
#  for ($i = 0; $i < @SampledFbgnAoA; $i++){
#    $k += @{$SampledFbgnAoA[$i]};
#  }
#  print "Random sites match $k total genes... ";
  
#Eliminate duplicate genes in SampledFbgnAoA (remove from the outlier region containing more genes)
  for ($i = 0; $i < @SampledFbgnAoA; $i++){
    for ($j = 0; $j < @{$SampledFbgnAoA[$i]}; $j++){
      $match = 0;
      for ($k = $i + 1; $k < @SampledFbgnAoA; $k++){
	for ($l = 0; $l < @{$SampledFbgnAoA[$k]}; $l++){
	  if ($SampledFbgnAoA[$i][$j] eq $SampledFbgnAoA[$k][$l]){
	    if (@{$SampledFbgnAoA[$i]} < @{$SampledFbgnAoA[$k]}){
	      @line = @{$SampledFbgnAoA[$k]};
	      splice @line, $l, 1;
	      @{$SampledFbgnAoA[$k]} = @line;
	      last;
	    }
	    else{
	      @line = @{$SampledFbgnAoA[$i]};
	      splice @line, $j, 1;
	      @{$SampledFbgnAoA[$i]} = @line;
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
#  $k = 0;
#  for ($i = 0; $i < @SampledFbgnAoA; $i++){
#    $k += @{$SampledFbgnAoA[$i]};
#  }
#  print "Matching $k unique genes... ";

#get GO categories from sampled regions
  for ($i = 0; $i < @SampledFbgnAoA; $i++){
    @RegionGOCats = ();
    for ($j = 0; $j < @{$SampledFbgnAoA[$i]}; $j++){
      $Fbgn = $SampledFbgnAoA[$i][$j];
      $Fbgn =~ s/FBgn//;
      $Fbgn += 0;
      for ($k = 0; $k < @{$GOAoA[$Fbgn]}; $k++){
	push @RegionGOCats, $GOAoA[$Fbgn][$k];
      }
    }
    @RegionGOCats  = sort{ lc($a) cmp lc($b) } @RegionGOCats;
    for ($j = 1; $j < @RegionGOCats; $j++){
      if ($RegionGOCats[$j] eq $RegionGOCats[$j-1]){
	splice @RegionGOCats, $j, 1;
	$j--;
      }
    }
    push @SampledGOList, @RegionGOCats;
  }
#  $k = @SampledGOList;
#  print "Associated with $k GO functions... ";
    
  for ($i = 0; $i < @AllGOUnique; $i++){
    push @SampledGOCounts, 0;
  }
  @SampledGOList = sort{ lc($a) cmp lc($b) } @SampledGOList;    
  $j = 0;
  for ($i = 0; $i < @SampledGOList; $i++){
    while ($SampledGOList[$i] ne $AllGOUnique[$j]){
      $j++;
    }
    $SampledGOCounts[$j]++;
  }

  
#Use a P-value array to add 1/reps each time the resampled data has more counts for a GO category than the empirical data did.
  for ($i = 0; $i < @OutlierGOCounts; $i++){
    if ($SampledGOCounts[$i] >= $OutlierGOCounts[$i]){
      $PValues[$i] += (1/$PermutationReps);
    }
  }
  print "Done with resampled data set $r\n";
}
close S;

#For each GO category with one or more outliers, look up the common gene names of the outliers
#@CdsFbgnAoA
#@CdsGeneAoA
my $GOCatOutliers = '';
my @GOCatOutlierList = ();
@OutlierFbgns = ();
my @OutlierGenes = ();
my @OutlierGOAoA = ();
for ($i = 0; $i < @OutlierFbgnAoA; $i++){
  for ($j = 0; $j < @{$OutlierFbgnAoA[$i]}; $j++){
    push @OutlierFbgns, $OutlierFbgnAoA[$i][$j];
  }
}
for ($i = 0; $i < @OutlierFbgns; $i++){
  $match = 0;
  for ($j = 0; $j < @CdsFbgnAoA; $j++){
    last if ($match == 1);
    for ($k = 0; $k < @{$CdsFbgnAoA[$j]}; $k++){
      if ($OutlierFbgns[$i] eq $CdsFbgnAoA[$j][$k]){
	push @OutlierGenes, $CdsGeneAoA[$j][$k];
	$Fbgn = $CdsFbgnAoA[$j][$k];
	$Fbgn =~ s/FBgn//;
	$Fbgn += 0;
	@line = @{$GOAoA[$Fbgn]};
	push @OutlierGOAoA, [ @line ];
	$match = 1;
	last;
      }
      if (($k == (@{$CdsFbgnAoA[$j]} - 1)) && ($j == (@CdsFbgnAoA - 1))){
	print "could not find $OutlierFbgns[$i] in CdsFbgnAoA\n";
      }
    }
  }	
}
$i = @OutlierFbgns;
$j = @OutlierGOAoA;
if ($i != $j){
  die "wrong number of entries in OutlierGOAoA ($j instead of $i)\n";
}
for ($i = 0; $i < @OutlierGOCounts; $i++){
  $GOCatOutliers = '';
  if ( ($OutlierGOCounts[$i] > 0) && ($i < (@OutlierGOCounts - 1))){
    for ($j = 0; $j < @OutlierGOAoA; $j++){
      for ($k = 0; $k < @{$OutlierGOAoA[$j]}; $k++){
	if ($AllGOUnique[$i] eq $OutlierGOAoA[$j][$k]){
	  next if ( defined($GOCatOutliers) && ($GOCatOutliers =~ m/$OutlierGenes[$j]/));
	  if ( defined($GOCatOutliers) && (length($GOCatOutliers)) > 0 ){
	    $GOCatOutliers = $GOCatOutliers . ',' . $OutlierGenes[$j];
	  }
	  else{
	    $GOCatOutliers = $OutlierGenes[$j];
	  }
	}
      }
    } 
    push @GOCatOutlierList, $GOCatOutliers;
  }
  else{
    $_ = '';
    push @GOCatOutlierList, $_;
  }
}

#Look up descriptions for each GO category
my @ontologies = ();
my @descriptions = ();
my @GODescAoA = ();
my $desc = '';
open D, "<$GOCatDescFile" or die "can not open $GOCatDescFile\n";
while (<D>){
  chomp;
  last if m/^$/;
  @line = split;
  push @GODescAoA, [ @line ];
}
close D;
for ($i = 0; $i < @AllGOUnique; $i++){
  for ($j = 0; $j < @GODescAoA; $j++){
    if ($AllGOUnique[$i] eq $GODescAoA[$j][0]){
      push @ontologies, $GODescAoA[$j][1];
      $desc = $GODescAoA[$j][2];
      for ($k = 3; $k < @{$GODescAoA[$j]}; $k++){
	$desc = $desc . ' ' . $GODescAoA[$j][$k];
      }
      push @descriptions, $desc;
      last;
    }
    if ($j == (@GODescAoA - 1)){
      $_ = '';
      push @ontologies, $_;
      push @descriptions, $_;
    }
  }
}

#Put all output in one AoA, sort by P value, and send to output file
my @OutputAoA = ();
for ($i = 0; $i < @AllGOUnique; $i++){
  @line = ();
  push @line, $AllGOUnique[$i];
  push @line, $ontologies[$i];
  push @line, $descriptions[$i];
  push @line, $OutlierGOCounts[$i];
  push @line, $AllGOCounts[$i];
  push @line, $PValues[$i];
  push @line, $GOCatOutlierList[$i];
  push @OutputAoA, [ @line ];
}
#@OutputAoA = sort {$a->[5] cmp $b->[5]} @OutputAoA;
open O, ">$OutputFile";
for ($i = 0; $i < @OutputAoA; $i++){
  for ($j = 0; $j < @{$OutputAoA[$i]}; $j++){
    print O $OutputAoA[$i][$j];
    if ($j < (@{$OutputAoA[$i]} - 1)){
      print O "\t";
    }
    else{
      print O "\n";
    }
  }
}
close O;
