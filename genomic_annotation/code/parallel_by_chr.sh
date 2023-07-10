#!/bin/bash

### 1. set paths ###

## set parameter parsings for user-specified paths ##
usage="This script parallelly runs genomic reannotations.\n
usage: sh parallel_by_chr [-h/--help] [-w/--work_directory <path>] [-s/--script <path>]\n
[-b/--bedtools <path>] [-a/--input_anno <path>] [-g/--input_genome_seq <path>]\n
[-c/--input_condon_table <path>] [-d/--conda_path <path>] [-e/--conda_env <name>]"

while [ $# -gt 0 ]; do
  case "$1" in
    --work_directory*|-w*)
      if [[ "$1" != *=* ]]; then shift; fi # Value is next arg if no `=`
      wd="${1#*=}"
      ;;
    --script*|-s*)
      if [[ "$1" != *=* ]]; then shift; fi
      script="${1#*=}"
      ;;
    --bedtools*|-b*)
      if [[ "$1" != *=* ]]; then shift; fi
      bedtools="${1#*=}"
      ;;
    --input_anno*|-a*)
      if [[ "$1" != *=* ]]; then shift; fi
      in_anno="${1#*=}"
      ;;
    --input_genome_seq*|-g*)
      if [[ "$1" != *=* ]]; then shift; fi
      in_genome_seq="${1#*=}"
      ;;
    --input_condon_table*|-c*)
      if [[ "$1" != *=* ]]; then shift; fi
      in_condon_tab="${1#*=}"
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
      echo $usage # Flag argument
      exit 0
      ;;
    *)
      >&2 printf "Error: Invalid argument\n"
      exit 1
      ;;
  esac
  shift
done

if [[ -z $wd || -z $script || -z $bedtools || -z $in_anno || -z $in_genome_seq || -z $in_condon_tab || -z $conda_path || -z $conda_env ]]; then
    echo $usage
    echo "-w: $wd
-s $script
-b $bedtools
-a $in_anno
-g $in_genome_seq
-c $input_condon_table
-d $conda_path
-e $conda_env"
    exit 1
fi

## derived paths ##
log=$wd/reannotate.log
temp=$wd/temp

sep_dir=$wd/sep
out_coord_sep=$sep_dir/reannotate_{}.bed
out_seq_sep_for=$sep_dir/for_reannotate_{}.fa
out_seq_sep_rev=$sep_dir/rev_reannotate_{}.fa
out_report_sep=$sep_dir/reannotate_{}.report
out_coord=$wd/reannotate.bed
out_seq_for=$wd/for_reannotate.fa
out_seq_rev=$wd/rev_reannotate.fa
out_report=$wd/reannotate.report

mkdir -p $sep_dir
mkdir -p $temp

### 2. configure the conda enviroment ###
source $conda_path/bin/activate $conda_env

### 3. run steps of the pipeline ###

### script starts ###
echo start at time `date +%F'  '%H:%M`

## define a function to reannotate each chromosome separately
function anno() {
local in_anno=$1
local in_genome_seq=$2
local in_condon_tab=$3
local out_coord_sep=$4
local out_seq_sep_for=$5
local out_seq_sep_rev=$6
local out_report_sep=$7
local bedtools=$8
local script=$9
local log=${10}
local chr=${11}

python $script \
--region $chr \
--in_anno $in_anno \
--in_genome_seq $in_genome_seq \
--in_condon_tab $in_condon_tab \
--out_coord $out_coord_sep \
--out_seq_for $out_seq_sep_for \
--out_seq_rev $out_seq_sep_rev \
--out_report $out_report_sep \
--bedtools $bedtools

echo "Reannotation done for $chr" >> $log
}

## function to parallel by chromosome/contig
function para() {
zcat $in_anno \
| awk 'BEGIN{FS="\t"} $1 !~ /#/{print $1}' \
| sort | uniq \
| parallel --tmpdir $temp -j 0 \
"$@"
}

## 3.1. run reannotations

echo -e "### Reannotations starts ###\n" >> $log

export -f anno

para \
anno $in_anno $in_genome_seq $in_condon_tab $out_coord_sep $out_seq_sep_for $out_seq_sep_rev $out_report_sep $bedtools $script $log {}

echo -e "### Reannotations ends ###\n" >> $log

## 3.2. merge outputs
echo -e "### Merging starts ###\n" >> $log

# merge bed files
para \
cat $out_coord_sep >> $out_coord

# merge fasta files
para \
cat $out_seq_sep_for >> $out_seq_for

para \
cat $out_seq_sep_rev >> $out_seq_rev

# merge report files
function ak() {
awk 'BEGIN{FS="\t"} \
NR==2 {printf "%s", $2} \
NR>2 {printf "\t%s", $2} \
END{printf "\n"}' $1
}

export -f ak
para ak $out_report_sep \
| awk 'BEGIN{FS="\t"; total=0} \
NR == 1 {print} \
NR > 1 && NR == FNR {name[NR - 1] = $1; len = NR-1} \
NR != FNR {for (i = 1; i <= NF; i++) a[i] += $i; total += $i} \
END{for (i = 1; i <= len; i++) printf "%s\t%d\t%-4.2f%%\n", name[i], a[i], a[i]*4*100/total}' \
<(zcat $in_anno | awk 'BEGIN{FS="\t"} $1 !~ /#/ {print $1}'| head -1 | parallel cat $out_report_sep) - \
> $out_report

echo -e "### Merging ends ###\n" >> $log

### script ends ###
echo finish at time `date +%F'  '%H:%M`