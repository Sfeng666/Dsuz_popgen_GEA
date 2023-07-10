#!/bin/bash

function run_paratest(){
mincov=$1
min_count=$2
window_heter=$3
wd=/home/siyuan/jobs/suzukii_WGS/genetic_diff/Fst/syn_Fst/results/Fst
sd=/home/siyuan/jobs/suzukii_WGS/genetic_diff/Fst/scripts
script=$sd/filter_mpileup_by_cov_window.py
calc_script=$sd/calc_Fst_reynolds_ungapped_efs.py
dir_snp=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/call_snp_maxcov_99/mincov_$mincov\_mincount_$min_count
dir_syn_snp=/home/siyuan/jobs/suzukii_WGS/genetic_diff/Fst/syn_Fst/results/syn_SNP_data
dir_mpileup=/home/siyuan/jobs/suzukii_WGS/06mpileup_rsync
chr_assignment=/home/siyuan/jobs/suzukii_WGS/calc_maxcov/amb/assignment_cor_amb.txt
sample_size=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/calc_diversity_maxcov_99/test_normalize/sample_size.txt
maxcov=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/calc_maxcov_99/max_cov.txt
conda_path=/opt/miniconda3
conda_env=WGS_analysis

bash /home/siyuan/jobs/suzukii_WGS/genetic_diff/Fst/syn_Fst/scripts/calc_syn_Fst_reynolds.sh \
-w $wd \
-s $dir_snp \
-y $dir_syn_snp \
-m $dir_mpileup \
-c $chr_assignment \
-a $sample_size \
-g $window_heter \
-v $mincov \
-x $maxcov \
-t $min_count \
-l $calc_script \
-p $script \
-d $conda_path \
-e $conda_env \
> $wd/../run_calc_syn_Fst_reynolds.log 2>&1
}

run_paratest 12 1 100