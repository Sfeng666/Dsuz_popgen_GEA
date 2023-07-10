#!/bin/bash

### 1. set paths ###
## set parameter parsings for user-specified paths ##
usage="## this script calculates biallelic minor allele frequency.\n
usage: bash $0 [-h/--help] [-w/--work_directory <path>] [-s/--dir_snp <path>]\n
[-m/--dir_mpileup <path>] [-c/--chr_assignment <path>]\n
[-v/--mincov <path>] [-t/--min_count <path>]\n
[-l/--calc_script <path>] [-d/--conda_path <path>] [-e/--conda_env <name>]"

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
    --min_count*|-t*)
      if [[ "$1" != *=* ]]; then shift; fi
      min_count="${1#*=}"
      ;;
    --calc_script*|-l*)
      if [[ "$1" != *=* ]]; then shift; fi
      calc_script="${1#*=}"
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

if [[ -z $wd || -z $dir_snp || -z $dir_mpileup \
|| -z $mincov || -z $min_count || -z $calc_script \
|| -z $chr_assignment || -z $conda_path || -z $conda_env ]]; then
    echo -e $usage
    echo "-w $wd
-s $dir_snp
-m $dir_mpileup
-c $chr_assignment
-v $mincov
-t $min_count
-l $calc_script
-d $conda_path
-e $conda_env"
    exit 1
fi

## derived paths ##
log=$wd/calc_maf.log
maf_auto=$wd/maf_auto.txt
maf_X=$wd/maf_X.txt
temp=$wd/temp
snp_clean_contig=$dir_snp/{}_snp_clean.vcf

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

calc_maf_contig() {
local snp_clean_contig=$1
local calc_script=$2
local log=$3

python $calc_script \
--in_snp $snp_clean_contig \

echo "maf calculated from $snp_clean_contig" >> $log
}

## For each contig,

## 3.1. calculate the minor allele frequency for InDel-free SNPs in each sample, then output the values into a record file
# note that the distributions of MAF are based on sites within autosomes-linked and X-linked (X-linked and ambiguous) contigs separately.

export -f calc_maf_contig

# autosome-linked contigs
echo -e "### Autosome-linked contigs: Estimation of dxy starts ###\n" >> $log
awk 'BEGIN{FS="_"} NR==1{printf $1} NR>1{printf "\t%s", $1} END{printf "\n"}' $bam_list > $maf_auto

{ time awk '$3 == "auto-linked"{print $1}' $chr_assignment \
| parallel --tmpdir $temp -j 0 \
calc_maf_contig $snp_clean_contig $calc_script $log >> $maf_auto ;} 2>> $log

gzip $maf_auto
echo -e "### Autosome-linked contigs: Estimation of dxy ends ###\n" >> $log

# X-linked and ambiguous
echo -e "### X-linked contigs: Estimation of dxy starts ###\n" >> $log
awk 'BEGIN{FS="_"} NR==1{printf $1} NR>1{printf "\t%s", $1} END{printf "\n"}' $bam_list > $maf_X

{ time awk '$3 == "X-linked" || $3 == "ambiguous"{print $1}' $chr_assignment \
| parallel --tmpdir $temp -j 0 \
calc_maf_contig $snp_clean_contig $calc_script $log >> $maf_X ;} 2>> $log

gzip $maf_X
echo -e "### X-linked contigs: Estimation of dxy ends ###\n" >> $log

echo -e "### Whole job done ###\n" >> $log
rm -r $temp

### script finishes ###
echo finish at time `date +%F'  '%H:%M`