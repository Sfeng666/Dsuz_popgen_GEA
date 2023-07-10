#!/bin/bash

### 1. set paths ###
## set parameter parsings for user-specified paths ##
usage="## this script calculates genome-wide window-based pairwise Reynolds Fst.\n
usage: bash $0 [-h/--help] [-w/--work_directory <path>] [-s/--dir_snp <path>]\n
[-m/--dir_syn_snp <path>] [-c/--chr_assignment <path>] [-a/--sample_size <path>]\n
[-g/--window_heter <int>] [-v/--mincov <int>] [-x/--maxcov <path>] [-t/--min_count <int>]\n
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
    --dir_syn_snp*|-y*)
      if [[ "$1" != *=* ]]; then shift; fi
      dir_syn_snp="${1#*=}"
      ;;
    --dir_mpileup*|-m*)
      if [[ "$1" != *=* ]]; then shift; fi
      dir_mpileup="${1#*=}"
      ;;      
    --chr_assignment*|-c*)
      if [[ "$1" != *=* ]]; then shift; fi
      chr_assignment="${1#*=}"
      ;;
    --sample_size*|-a*)
      if [[ "$1" != *=* ]]; then shift; fi
      sample_size="${1#*=}"
      ;;
    --window_heter*|-g*)
      if [[ "$1" != *=* ]]; then shift; fi
      window_heter="${1#*=}"
      ;;
    --min_cov*|-v*)
      if [[ "$1" != *=* ]]; then shift; fi
      min_cov="${1#*=}"
      ;;
    --max_cov*|-x*)
      if [[ "$1" != *=* ]]; then shift; fi
      max_cov="${1#*=}"
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

if [[ -z $wd || -z $dir_snp || -z $dir_syn_snp || -z $dir_mpileup || -z $chr_assignment  \
|| -z $sample_size || -z $window_heter \
|| -z $max_cov || -z $min_cov || -z $min_count || -z $calc_script || -z $script \
|| -z $conda_path || -z $conda_env ]]; then
    echo -e $usage
    echo "-w $wd
-s $dir_snp
-y $dir_syn_snp
-m $dir_mpileup
-c $chr_assignment
-a $sample_size
-g $window_heter
-v $min_cov
-x $max_cov
-t $min_count
-l $calc_script
-p $script
-d $conda_path
-e $conda_env"
    exit 1
fi

## derived paths ##
log=$wd/calc_Fst.log
Fst_auto=$wd/Fst_auto.txt
Fst_X=$wd/Fst_X.txt
window_stats=$wd/window_stats.txt
temp=$wd/temps
snp_clean_contig=$dir_syn_snp/{1}_snp_clean_syn.vcf
indel_positions_contig=$dir_snp/{1}_inDel-positions_20.txt
mpileup_contig=$dir_syn_snp/{1}_syn_mpileup.gz

bam_list=$dir_mpileup/bam_list.txt

mkdir -p $wd
mkdir -p $temp

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" 'Chr' 'Contig' 'Window' 'Start' 'End' 'Length' 'Number of filtered sites' 'Number of SNPs' > $window_stats

names=""
while read -r line
do
    sample=$(basename $line _sort_merge_dedup_indel.bam)
    if [ $names ]
    then
        names+=",${sample}"
    else
        names+=$sample
    fi
done < $bam_list

# ### 2. configure the conda enviroment ###
# source $conda_path/bin/activate $conda_env

### 3. run steps of the pipeline ###
### job starts ###
echo start at time `date +%F'  '%H:%M`

echo -e "### Whole job starts ###\n" > $log

## define a function to calculate the numerator and denominator of window-based pairwise Fst
calc_Fst_contig() {
local contig=$1
local chr=$2
local snp_clean_contig=$3
local indel_positions_contig=$4
local mpileup_contig=$5

python $calc_script \
--in_snp $snp_clean_contig \
--in_mpileup $mpileup_contig \
--in_indel_pos $indel_positions_contig \
--in_sample_size $sample_size \
--out_window_stats $window_stats \
--script $script \
--window_heter $window_heter \
--min_count $min_count \
--min_cov $min_cov \
--max_cov $max_cov \
--contig $contig \
--chr $chr \
--names $names \
--temp $temp

echo "Fst calculated from $snp_clean_contig" >> $log
}

## For each contig,

## 3.1. calculate the accumulative weighted pairwise window-Fst for InDel-free SNPs, then devided by the genome-wide total number of filtered sites.
# note that the calculations are based on sites within autosomes-linked and X-linked (X-linked and ambiguous) contigs separately.

export -f calc_Fst_contig
export sample_size
export window_stats
export calc_script
export script
export window_heter
export min_count
export min_cov
export max_cov
export names
export temp
export log

# autosome-linked contigs
echo -e "### Autosome-linked contigs: Estimation of Fst starts ###\n" >> $log

{ time awk 'BEGIN{OFS="\t"} $3 == "auto-linked"{print $1, "auto"}' $chr_assignment \
| parallel --tmpdir $temp -j 0 --colsep '\t' \
calc_Fst_contig {1} {2} $snp_clean_contig $indel_positions_contig $mpileup_contig \
| awk 'BEGIN{FS = OFS = "\t"} \
NR==1{for (i = 2; i <= NF; i++) {samples[i] = $i}; num_sites += $1; print} \
NR>1 && $1 ~/[0-9]+/{num_sites += $1} \
NR>1 && $1 !~/[0-9]+/{sample=$1; for (i=2; i<=NF; i++) Fst[sample, i] += $i} \
END{for (i=2; i<=NF; i++) {printf "%s\t", samples[i]; for (j=2; j<=NF-1; j++) {printf "%f\t", Fst[samples[i], j]/num_sites}; printf "%f\n", Fst[samples[i], j]/num_sites}}' \
1> $Fst_auto  ;} 2>> $log

echo -e "### Autosome-linked contigs: Estimation of Fst ends ###\n" >> $log

# X-linked and ambiguous
echo -e "### X-linked contigs: Estimation of Fst starts ###\n" >> $log

{ time awk 'BEGIN{OFS="\t"} $3 == "X-linked" || $3 == "ambiguous"{print $1, "X"}' $chr_assignment \
| parallel --tmpdir $temp -j 0 --colsep '\t' \
calc_Fst_contig {1} {2} $snp_clean_contig $indel_positions_contig $mpileup_contig \
| awk 'BEGIN{FS = OFS = "\t"} \
NR==1{for (i = 2; i <= NF; i++) {samples[i] = $i}; num_sites += $1; print} \
NR>1 && $1 ~/[0-9]+/{num_sites += $1} \
NR>1 && $1 !~/[0-9]+/{sample=$1; for (i=2; i<=NF; i++) Fst[sample, i] += $i} \
END{for (i=2; i<=NF; i++) {printf "%s\t", samples[i]; for (j=2; j<=NF-1; j++) {printf "%f\t", Fst[samples[i], j]/num_sites}; printf "%f\n", Fst[samples[i], j]/num_sites}}' \
1> $Fst_X  ;} 2>> $log

echo -e "### X-linked contigs: Estimation of Fst ends ###\n" >> $log

echo -e "### Whole job done ###\n" >> $log
rm -r $temp

### script finishes ###
echo finish at time `date +%F'  '%H:%M`