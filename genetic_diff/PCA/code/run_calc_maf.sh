#!/bin/bash

function run_paratest(){
mincov=$1
min_count=$2
wd=/home/siyuan/jobs/suzukii_WGS/genetic_diff/PCA/results/mincov_$mincov\_mincount_$min_count
dir_snp=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/call_snp_maxcov_99/mincov_$mincov\_mincount_$min_count
dir_mpileup=/home/siyuan/jobs/suzukii_WGS/06mpileup_rsync
calc_script=/home/siyuan/jobs/suzukii_WGS/genetic_diff/PCA/corrected_sortallele/calc_maf_biallelic.py
chr_assignment=/home/siyuan/jobs/suzukii_WGS/calc_maxcov/amb/assignment_cor_amb.txt
conda_path=/opt/miniconda3
conda_env=WGS_analysis

bash calc_maf.sh \
-w $wd \
-s $dir_snp \
-m $dir_mpileup \
-c $chr_assignment \
-v $mincov \
-t $min_count \
-l $calc_script \
-d $conda_path \
-e $conda_env \
> run_mincov_$mincov\_mincount_$min_count\.log 2>&1
}

run_paratest 12 1
# run_paratest 12 10
# run_paratest 12 5
# run_paratest 12 8
# run_paratest 12 3
# run_paratest 12 20
# run_paratest 12 50
# run_paratest 20 10