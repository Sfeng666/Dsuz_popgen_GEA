#!/bin/bash

### 1. set paths ###
## set parameter parsings for user-specified paths ##
usage="## this script builds maximum-likelihood tree using TreeMix.\n
usage: bash $0 [-h/--help] [-w/--work_directory <path>] [-s/--dir_snp <path>]\n
[-c/--chr_assignment <path>] [-a/--sample_size <path>] [-l/--window_len <int>]\n
[-p/--script <path>] [-d/--conda_path <path>] [-e/--conda_env <name>]"

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
    --chr_assignment*|-c*)
      if [[ "$1" != *=* ]]; then shift; fi
      chr_assignment="${1#*=}"
      ;;
    --sample_size*|-a*)
      if [[ "$1" != *=* ]]; then shift; fi
      sample_size="${1#*=}"
      ;;
    --window_len*|-l*)
      if [[ "$1" != *=* ]]; then shift; fi
      window_len="${1#*=}"
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

if [[ -z $wd || -z $dir_snp || -z $chr_assignment \
|| -z $sample_size || -z $window_len || -z $script \
|| -z $conda_path || -z $conda_env ]]; then
    echo -e $usage
    echo "-w $wd
-s $dir_snp
-c $chr_assignment
-a $sample_size
-l $window_len
-p $script
-d $conda_path
-e $conda_env"
    exit 1
fi

## derived paths ##
log=$wd/calc_MLtree.log
in_auto=$wd/in_auto.txt
in_X=$wd/in_X.txt
in_any=$wd/in_{2}.txt
snpnum_auto=$wd/snpnum_auto.txt
snpnum_X=$wd/snpnum_X.txt
snpnum_any=$wd/snpnum_any.txt
winlen_contig=$wd/winlen_contig.txt
temp=$wd/tempsq
snp_clean_contig=$dir_snp/{1}_snp_clean.vcf

mkdir -p $wd
mkdir -p $temp

awk 'BEGIN{FS = "\t"; OFS = " "} NR == 1{print $0}' $sample_size > $in_auto 
awk 'BEGIN{FS = "\t"; OFS = " "} NR == 1{print $0}' $sample_size > $in_X 

### 2. configure the conda enviroment ###
source $conda_path/bin/activate $conda_env #somehow if you do not activate conda, a perl library error will occur at the parallel step

### 3. run steps of the pipeline ###
### job starts ###
echo start at time `date +%F'  '%H:%M`

echo -e "### Whole job starts ###\n" > $log

## define a function to generate the input table of SNP allele frequency for TreeMix, and generate window length for comparison
prep_treemix_input() {
local snp_clean_contig=$1
local chr=$2
local out_window_len=$3
local out_SNP_num=$4

python $script \
--in_snp $snp_clean_contig \
--in_sample_size $sample_size \
--out_window_len $out_window_len \
--out_SNP_num $out_SNP_num \
--chr $chr \
--temp $temp

echo "Allele frequency extracted from $snp_clean_contig" >> $log
}

## For each contig,

## 3.1. write the input table of SNP allele frequency for TreeMix, and generate window length for comparison.
# note that the extractions are based on sites within autosomes-linked and X-linked (X-linked and ambiguous) contigs separately.

export -f prep_treemix_input
export sample_size
export script
export in_auto
export in_X
export temp
export log

# 3.1.1. output SNP allele frequency and snp number of each contig to separate auto- and X-linked TreeMix input files
echo -e "### Auto-linked contigs: Prep of TreeMix input starts ###\n" >> $log

{ time awk 'BEGIN{OFS="\t"} $3 == "auto-linked"{print $1, "auto"}' $chr_assignment \
| parallel --tmpdir $temp -k -j 0 --colsep '\t' \
prep_treemix_input $snp_clean_contig {2} $winlen_contig $snpnum_auto >> $in_auto ;} 2>> $log

## for test only
# dir_snp=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/call_snp_maxcov_99/mincov_12_mincount_10
# chr_assignment=/home/siyuan/jobs/suzukii_WGS/calc_maxcov/amb/assignment_cor_amb.txt
# awk 'BEGIN{OFS="\t"} $3 == "auto-linked"{print $1, "auto"}' $chr_assignment \
# | parallel --keep-order -j 0 --colsep '\t' \
# tail -2 $dir_snp/{1}_snp_clean.vcf >> test

for contig in $(awk 'BEGIN{OFS="\t"} $3 == "auto-linked"{print $1}' $chr_assignment)  # keep the contig order of SNP number the same as the SNP array
do
  grep -v '#' $dir_snp/$contig\_snp_clean.vcf | wc -l >> $snpnum_auto
done

echo -e "### Auto-linked contigs: Prep of TreeMix input ends ###\n" >> $log

echo -e "### X-linked contigs: Prep of TreeMix input starts ###\n" >> $log

{ time awk 'BEGIN{OFS="\t"} $3 == "X-linked" || $3 == "ambiguous"{print $1, "X"}' $chr_assignment \
| parallel --tmpdir $temp -k -j 0 --colsep '\t' \
prep_treemix_input $snp_clean_contig {2} $winlen_contig $snpnum_X >> $in_X ;} 2>> $log

for contig in $(awk 'BEGIN{OFS="\t"} $3 == "X-linked" || $3 == "ambiguous"{print $1}' $chr_assignment) # keep the contig order of SNP number the same as the SNP array
do
  grep -v '#' $dir_snp/$contig\_snp_clean.vcf | wc -l >> $snpnum_X
done

echo -e "### X-linked contigs: Prep of TreeMix input ends ###\n" >> $log

# 3.1.2. calculate K across all contigs
{ time K=$(awk -v window_len=$window_len \
'BEGIN{FS = OFS = "\t"; min_diff = window_len}
{K = $1; for (i = 2; i <= NF; i++) {all_len[K][NR, i] = $i; win_num[K] += 1}}
END{for (K in all_len) {asort(all_len[K])
median_len = (win_num[K] %2 == 0) ? (all_len[K][win_num[K]/2] + all_len[K][(win_num[K]/2) + 1])/2 : all_len[K][int(win_num[K]/2) + 1]
diff = sqrt((median_len - window_len)^2); if (diff < min_diff) {min_diff = diff; min_K = K}}; print min_K}' $winlen_contig)  ;} 2>> $log

echo "-K for all contigs = $K" >> $log
echo -e "### All contigs: Prep of TreeMix input ends ###\n" >> $log

# 3.1.3. extract SNPs based on the determeined window size, and discard 3' SNPs as the remainder

## annnotated below is a problematic SNP extraction - cannot correctly deal with the case where the number of SNPs within a contig could be whole-divided by K, and will lead to a shift in following window partition 
# echo -e "### Autosome-linked contigs: extraction of windowed SNPs starts ###\n" >> $log
# awk -v K=$K 'BEGIN{FS = OFS = " "; contig = contig_crt = pos = 1} 
# ARGIND == 1 && $2 != 0{snp_num[contig] = $2; contig += 1}
# ARGIND == 2 && FNR == 1 {print}
# ARGIND == 2 && FNR > 1{if (pos <= snp_num[contig_crt] - snp_num[contig_crt]%K) {print; pos += 1; if (snp_num[contig_crt]%K == 0) {contig_crt += 1; pos = 1}}
# else if (pos == snp_num[contig_crt]) {contig_crt += 1; pos = 1} 
# else pos += 1}' $snpnum_auto $in_auto \
# | gzip > $in_auto\.gz
# echo -e "### Autosome-linked contigs: extraction of windowed SNPs ends ###\n" >> $log

echo -e "### Autosome-linked contigs: extraction of windowed SNPs starts ###\n" >> $log
awk -v K=$K 'BEGIN{FS = OFS = " "; contig = contig_crt = pos = 1} 
ARGIND == 1 && $1 != 0{snp_num[contig] = $1; contig += 1}
ARGIND == 2 && FNR == 1 {print}
ARGIND == 2 && FNR > 1{if (snp_num[contig_crt]%K == 0) {if (pos < snp_num[contig_crt]) {print; pos += 1} else if (pos == snp_num[contig_crt]) {print; contig_crt += 1; pos = 1}}
else if (snp_num[contig_crt]%K > 0) {if (pos <= snp_num[contig_crt] - snp_num[contig_crt]%K) {print; pos += 1} else if (pos > snp_num[contig_crt] - snp_num[contig_crt]%K) {if (pos < snp_num[contig_crt]) {pos += 1} else if (pos == snp_num[contig_crt]) {contig_crt += 1; pos = 1}}}}' \
$snpnum_auto $in_auto \
| gzip > $in_auto\.gz
echo -e "### Autosome-linked contigs: extraction of windowed SNPs ends ###\n" >> $log

## annnotated below is a problematic SNP extraction - cannot correctly deal with the case where the number of SNPs within a contig could be whole-divided by K, and will lead to a shift in following window partition 
# echo -e "### X-linked contigs: extraction of windowed SNPs starts ###\n" >> $log
# awk -v K=$K 'BEGIN{FS = OFS = " "; contig = contig_crt = pos = 1}
# ARGIND == 1 && $2 != 0{snp_num[contig] = $2; contig += 1}
# ARGIND == 2 && FNR == 1 {print}
# ARGIND == 2 && FNR > 1{if (pos <= snp_num[contig_crt] - snp_num[contig_crt]%K) {print; pos += 1; if (snp_num[contig_crt]%K == 0) {contig_crt += 1; pos = 1}}
# else if (pos == snp_num[contig_crt]) {contig_crt += 1; pos = 1}
# else pos += 1}' $snpnum_X $in_X \
# | gzip > $in_X\.gz
# echo -e "### X-linked contigs: extraction of windowed SNPs ends ###\n" >> $log

echo -e "### X-linked contigs: extraction of windowed SNPs starts ###\n" >> $log
awk -v K=$K 'BEGIN{FS = OFS = " "; contig = contig_crt = pos = 1} 
ARGIND == 1 && $1 != 0{snp_num[contig] = $1; contig += 1}
ARGIND == 2 && FNR == 1 {print}
ARGIND == 2 && FNR > 1{if (snp_num[contig_crt]%K == 0) {if (pos < snp_num[contig_crt]) {print; pos += 1} else if (pos == snp_num[contig_crt]) {print; contig_crt += 1; pos = 1}}
else if (snp_num[contig_crt]%K > 0) {if (pos <= snp_num[contig_crt] - snp_num[contig_crt]%K) {print; pos += 1} else if (pos > snp_num[contig_crt] - snp_num[contig_crt]%K) {if (pos < snp_num[contig_crt]) {pos += 1} else if (pos == snp_num[contig_crt]) {contig_crt += 1; pos = 1}}}}' \
$snpnum_X $in_X \
| gzip > $in_X\.gz
echo -e "### X-linked contigs: extraction of windowed SNPs ends ###\n" >> $log

## 3.2. run TreeMix with a set of different number of migration events (10)
run_treemix_mig() {
local chr=$1
local out_treemix_in=$2
local m=$3

echo TreeMix at $chr with $m migration events starts at time `date +%F'  '%H:%M`  >> $log
treemix -i $out_treemix_in\.gz -k $K -root CN-Nin -m $3 -o $wd/$chr\_m_$m >> $log 2>&1
echo TreeMix at $chr with $m migration events ends at time `date +%F'  '%H:%M`  >> $log
}

export -f run_treemix_mig
export K
export wd

echo -e "### X-linked contigs: Running of TreeMix starts ###\n" >> $log
treemix -i $in_X\.gz -k $K -root CN-Nin -o $wd/X_m_0
seq 0 20 \
| parallel --memfree 10G --tmpdir $temp -j 0 \
run_treemix_mig X $in_X {} &

echo -e "### X-linked contigs: Running of TreeMix ends ###\n" >> $log

echo -e "### Auto-linked contigs: Running of TreeMix starts ###\n" >> $log
treemix -i $in_auto\.gz -k $K -root CN-Nin -o $wd/auto_m_0
seq 0 20 \
| parallel --memfree 10G --tmpdir $temp -j 0 \
run_treemix_mig auto $in_auto {}

echo -e "### Auto-linked contigs: Running of TreeMix ends ###\n" >> $log

echo -e "### Whole job done ###\n" >> $log
rm -r $temp

### script finishes ###
echo finish at time `date +%F'  '%H:%M`

