#!/bin/bash

### 1. set paths ###
## set parameter parsings for user-specified paths ##
usage="## this script calculates pairwise genetic distance dxy.\n
usage: bash $0 [-h/--help] [-w/--work_directory <path>] [-s/--dir_snp <path>]\n
[-m/--dir_mpileup <path>] [-c/--chr_assignment <path>]\n
[-v/--mincov <path>] [-x/--maxcov <path>] [-t/--min_count <path>]\n
[-l/--calc_script <path>] [-p/--script <path>] [-d/--conda_path <path>] [-e/--conda_env <name>]"

while [ $# -gt 0 ]; do
  case "$1" in
    --work_directory*|-w*)
      if [[ "$1" != *=* ]]; then shift; fi # Value is next arg if no `=`
      wd="${1#*=}"
      ;;
    --dir_snp*|-s*)
      if [[ "$1" != *=* ]]; then shift; fi
      dir_snp="${1#*=}"
      ;;
    --dir_mpileup*|-m*)
      if [[ "$1" != *=* ]]; then shift; fi
      dir_mpileup="${1#*=}"
      ;;
    --chr_assignment*|-c*)
      if [[ "$1" != *=* ]]; then shift; fi
      chr_assignment="${1#*=}"
      ;;
    --mincov*|-v*)
      if [[ "$1" != *=* ]]; then shift; fi
      mincov="${1#*=}"
      ;;
    --maxcov*|-x*)
      if [[ "$1" != *=* ]]; then shift; fi
      maxcov="${1#*=}"
      ;;
    --min_count*|-t*)
      if [[ "$1" != *=* ]]; then shift; fi
      min_count="${1#*=}"
      ;;
    --calc_script*|-l*)
      if [[ "$1" != *=* ]]; then shift; fi
      calc_script="${1#*=}"
      ;;
    --script*|-p*)
      if [[ "$1" != *=* ]]; then shift; fi
      script="${1#*=}"
      ;;
    --conda_path*|-d*)
      if [[ "$1" != *=* ]]; then shift; fi
      conda_path="${1#*=}"
      ;;
    --conda_env*|-e*)
      if [[ "$1" != *=* ]]; then shift; fi
      conda_env="${1#*=}"
      ;;      
    --help|-h)
      echo -e $usage # Flag argument
      exit 0
      ;;
    *)
      >&2 printf "Error: Invalid argument\n"
      exit 1
      ;;
  esac
  shift
done

if [[ -z $wd || -z $dir_snp || -z $dir_mpileup || -z $maxcov \
|| -z $mincov || -z $min_count || -z $calc_script || -z $script \
|| -z $chr_assignment || -z $conda_path || -z $conda_env ]]; then
    echo -e $usage
    echo "-w $wd
-s $dir_snp
-m $dir_mpileup
-c $chr_assignment
-v $mincov
-x $maxcov
-t $min_count
-l $calc_script
-p $script
-d $conda_path
-e $conda_env"
    exit 1
fi

## derived paths ##
log=$wd/calc_dxy.log
distance_auto=$wd/genetic_distance_auto.txt
distance_X=$wd/genetic_distance_X.txt
temp=$wd/temp
snp_clean_contig=$dir_snp/{}_snp_clean.vcf
indel_positions_contig=$dir_snp/{1}_inDel-positions_20.txt
mpileup_contig=$dir_mpileup/{1}_mpileup.gz
min_freq=0
miss_frac=0.001

bam_list=$dir_mpileup/bam_list.txt

mkdir -p $wd
mkdir -p $temp

names=""
num_sample=0
while read -r line
do
    sample=$(basename $line _sort_merge_dedup_indel.bam)
    if [ $names ]
    then
        names+=",${sample}"
    else
        names+=$sample
    fi
    ((num_sample++))
done < $bam_list

### 2. configure the conda enviroment ###
source $conda_path/bin/activate $conda_env

### 3. run steps of the pipeline ###
### job starts ###
echo start at time `date +%F'  '%H:%M`

echo -e "### Whole job starts ###\n" > $log

## define a function to calculate total heterozygosity for InDel-free SNPs within each contig for each sample

calc_dxy_contig() {
local snp_clean_contig=$1
local names=$2
local calc_script=$3
local log=$4

python $calc_script \
--in_snp $snp_clean_contig \
--names $names

echo "dxy calculated from $snp_clean_contig" >> $log
}

## define a function to count the number of InDel-free and coverage-filtered sites
count_sites() {
local indel_positions_contig=$1
local mpileup_contig=$2
local min_cov=$3
local max_cov=$4
local min_count=$5
local min_freq=$6
local miss_frac=$7
local names=$8
local script=$9
local log=${10}
local temp=${11}


awk 'BEGIN{FS="\t"; sites=0} \
NR==FNR {D[$1$2]++; next} \
!($1$2 in D) && NR>FNR {print} ' \
<(cat $indel_positions_contig) \
<(zcat $mpileup_contig) \
| parallel \
--tmpdir $temp \
-k \
--pipe \
-j 100% \
--no-notice \
--cat python2.7 $script \
--mpileup  {} \
--min-cov $min_cov \
--max-cov $max_cov \
--min-count $min_count \
--min-freq $min_freq \
--miss-frac $miss_frac \
--names $names

echo "InDel-free and coverage-filtered sites counted from $mpileup_contig" >> $log
}

## For each contig,
## 3.1. count the number of InDel-free filtered sites, then output the values into a record variable

# autosome-linked contigs
echo -e "### Autosome-linked contigs: Counting of InDel-free filtered sites starts ###\n" >> $log

export -f count_sites
{ time num_sites_auto=$(awk 'BEGIN{OFS="\t"} $4 == "auto-linked"{print $1, $2}' <(join -1 1 -2 1 <(sort -k 1 $maxcov) <(sort -k 1 $chr_assignment)) \
| parallel --tmpdir $temp -j 0 --colsep '\t' \
count_sites $indel_positions_contig $mpileup_contig $mincov $maxcov $min_count $min_freq $miss_frac $names $script $log $temp \
| awk -v num_sample=$num_sample 'BEGIN{sites=0} \
{sites+=$1} \
END{print sites}') ;} 2>> $log

echo "Autosome-linked contigs: number of InDel-free filtered sites = $num_sites_auto"  >> $log
echo -e "### Autosome-linked contigs: Counting of InDel-free filtered sites ends ###\n" >> $log

# X-linked and ambiguous
echo -e "### X-linked contigs: Counting of InDel-free filtered sites starts ###\n" >> $log

export -f count_sites
{ time num_sites_X=$(awk 'BEGIN{OFS="\t"} $4 == "X-linked" || $4 == "ambiguous"{print $1, $2}' <(join -1 1 -2 1 <(sort -k 1 $maxcov) <(sort -k 1 $chr_assignment)) \
| parallel --tmpdir $temp -j 0 --colsep '\t' \
count_sites $indel_positions_contig $mpileup_contig $mincov $maxcov $min_count $min_freq $miss_frac $names $script $log $temp \
| awk -v num_sample=$num_sample 'BEGIN{sites=0} \
{sites+=$1} \
END{print sites}') ;} 2>> $log

echo "X-linked contigs: number of InDel-free filtered sites = $num_sites_X"  >> $log
echo -e "### X-linked contigs: Counting of InDel-free filtered sites ends ###\n" >> $log

# # test only - test if part 3.2 works
# num_sites_auto=119650000
# num_sites_X=19865000

## 3.2. calculate the accumulative value of heterozygosity for InDel-free SNPs in each sample, then output the values into a record file.
# note that the calculations are based on sites within autosomes-linked and X-linked (X-linked and ambiguous) contigs separately.

export -f calc_dxy_contig

# autosome-linked contigs
echo -e "### Autosome-linked contigs: Estimation of dxy starts ###\n" >> $log

{ time awk '$4 == "auto-linked"{print $1}' <(join -1 1 -2 1 <(sort -k 1 $maxcov) <(sort -k 1 $chr_assignment)) \
| parallel --tmpdir $temp -j 0 \
calc_dxy_contig $snp_clean_contig $names $calc_script $log \
| awk 'BEGIN{FS=OFS="\t"} \
NR==1{num_sites=$1} \
NR==2{for (i=2; i<=NF; i++) {samples[i]=$i}; print} \
NR>2 && $0 !~/Populations/{sample=$1; for (i=2; i<=NF; i++) dxy[sample, i] += $i} \
END{for (i=2; i<=NF; i++) {printf "%s\t", samples[i]; for (j=2; j<=NF-1; j++) {printf "%f\t", dxy[samples[i], j]/num_sites}; printf "%f\n", dxy[samples[i], j]/num_sites}}' \
<(echo $num_sites_auto) - \
1> $distance_auto  ;} 2>> $log

echo -e "### Autosome-linked contigs: Estimation of dxy ends ###\n" >> $log

# X-linked and ambiguous
echo -e "### X-linked contigs: Estimation of dxy starts ###\n" >> $log

{ time awk '$4 == "X-linked" || $4 == "ambiguous"{print $1}' <(join -1 1 -2 1 <(sort -k 1 $maxcov) <(sort -k 1 $chr_assignment)) \
| parallel --tmpdir $temp -j 0 \
calc_dxy_contig $snp_clean_contig $names $calc_script $log \
| awk 'BEGIN{FS=OFS="\t"} \
NR==1{num_sites=$1} \
NR==2{for (i=2; i<=NF; i++) {samples[i]=$i}; print} \
NR>2 && $0 !~/Populations/{sample=$1; for (i=2; i<=NF; i++) dxy[sample, i] += $i} \
END{for (i=2; i<=NF; i++) {printf "%s\t", samples[i]; for (j=2; j<=NF-1; j++) {printf "%f\t", dxy[samples[i], j]/num_sites}; printf "%f\n", dxy[samples[i], j]/num_sites}}' \
<(echo $num_sites_X) - \
1> $distance_X  ;} 2>> $log

echo -e "### X-linked contigs: Estimation of dxy ends ###\n" >> $log

echo -e "### Whole job done ###\n" >> $log
rm -r $temp

### script finishes ###
echo finish at time `date +%F'  '%H:%M`