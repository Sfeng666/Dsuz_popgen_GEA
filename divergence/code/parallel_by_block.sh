#!/bin/bash

### 1. set paths ###
## set parameter parsings for user-specified paths ##
usage="## Parallely extract SNPs with known predicted ancestral alleles.\n
usage: sh $0 [-h/--help] [-w/--work_directory <path>] [-a/--in_anno_seq <path>]\n
[-m/--in_map <path>] [-s/--script <path>] [-r/--in_assembly_rep <path>]\n
[-t/--in_mfa_tgt <path>] [-n/--in_mfa_ans <path>]\n
[-b/--bedtools <path>] [-d/--conda_path <path>] [-e/--conda_env <name>]"

while [ $# -gt 0 ]; do
  case "$1" in
    --work_directory*|-w*)
      if [[ "$1" != *=* ]]; then shift; fi # Value is next arg if no `=`
      wd="${1#*=}"
      ;;
    --input_anno_seq*|-a*)
      if [[ "$1" != *=* ]]; then shift; fi
      in_anno_seq="${1#*=}"
      ;;
    --in_map*|-m*)
      if [[ "$1" != *=* ]]; then shift; fi
      in_map="${1#*=}"
      ;;
    --script*|-s*)
      if [[ "$1" != *=* ]]; then shift; fi
      script="${1#*=}"
      ;;
    --in_assembly_rep*|-r*)
      if [[ "$1" != *=* ]]; then shift; fi
      in_assembly_rep="${1#*=}"
      ;;
    --in_mfa_tgt*|-t*)
      if [[ "$1" != *=* ]]; then shift; fi
      in_mfa_tgt="${1#*=}"
      ;;
    --in_mfa_ans*|-n*)
      if [[ "$1" != *=* ]]; then shift; fi
      in_mfa_ans="${1#*=}"
      ;;
    --bedtools*|-b*)
      if [[ "$1" != *=* ]]; then shift; fi
      bedtools="${1#*=}"
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

if [[ -z $wd || -z $in_anno_seq || -z $in_map || -z $script || -z $in_assembly_rep \
|| -z $in_mfa_tgt || -z $in_mfa_ans || -z $bedtools || -z $conda_path || -z $conda_env ]]; then
    echo -e $usage
    echo "-w: $wd
-a $in_anno_seq
-m $in_map
-s $script
-r $in_assembly_rep
-t $in_mfa_tgt
-n $in_mfa_ans
-b $bedtools
-d $conda_path
-e $conda_env"
    exit 1
fi

## derived paths ##
in_mfa_tgt=$in_mfa_tgt\_{1}.mfa
in_mfa_ans=$in_mfa_ans\_{1}.mfa
out_divergence=$wd/divergence.txt
log=$wd/calc_divergence.log
temp=$wd/temp

mkdir -p $temp

### 2. configure the conda enviroment ###
source $conda_path/bin/activate $conda_env

### 3. run steps of the pipeline ###

### script starts ###
echo start at time `date +%F'  '%H:%M`

echo -e "### Whole job starts ###\n" > $log

## define a function to extract SNPs
function calc_div() {
local blk=$1
local in_assembly_rep=$2
local in_anno_seq=$3
local in_map=$4
local in_mfa_tgt=$5
local in_mfa_ans=$6
local temp=$7
local bedtools=$8
local script=$9
local chr=${10}
local log=${11}

python $script \
--blk $blk \
--in_assembly_rep $in_assembly_rep \
--in_anno_seq $in_anno_seq \
--in_map $in_map \
--in_mfa_tgt $in_mfa_tgt \
--in_mfa_ans $in_mfa_ans \
--temp $temp \
--bedtools $bedtools \
2>> $log

echo "Divergence calculated for $chr $blk" >> $log
}

## For each block,
## 3.1. extract SNPs with known predicted ancestral alleles.
echo -e "### Calculation of divergence starts ###\n" >> $log

export -f calc_div
awk 'BEGIN{FS = OFS = "\t"} \
$2 != "NA" && $6 != "NA"{print $1, $2}' $in_map \
| parallel --no-notice --tmpdir $temp -j 0 --colsep '\t' \
calc_div {1} $in_assembly_rep $in_anno_seq $in_map $in_mfa_tgt $in_mfa_ans $temp $bedtools $script {2} $log \
| awk 'BEGIN{FS = "\t"} \
NR == 1 {print} \
int(NR/2) == NR/2 {i = 1; while (i <= NF/2) {total[i] += $(2*i); div[i] += $(2*i - 1); i ++}} \
END{for (i=1;i<=NF/2-1;i+=1) printf "%f\t", div[i]/total[i]; printf "%f\n", div[i]/total[i]}' \
> $out_divergence

rm -r $temp
echo -e "### Calculation of divergence  ends ###\n" >> $log

echo -e "### Whole job done ###\n" >> $log

### script finishes ###
echo finish at time `date +%F'  '%H:%M`
