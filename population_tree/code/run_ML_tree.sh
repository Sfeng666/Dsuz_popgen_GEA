#!/bin/bash

function run_paratest(){
mincov=$1
min_count=$2
window_len=$3
wd=/home/siyuan/jobs/suzukii_WGS/genetic_diff/ML_tree/results/mincount_10_no_se_auto_X_over10/winlen_$window_len\_mincov_$mincov\_mincount_$min_count
sd=/home/siyuan/jobs/suzukii_WGS/genetic_diff/ML_tree/scripts
script=$sd/prep_input_forTreeMix_efs.py
dir_snp=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/call_snp_maxcov_99/mincov_$mincov\_mincount_$min_count
chr_assignment=/home/siyuan/jobs/suzukii_WGS/calc_maxcov/amb/assignment_cor_amb.txt
sample_size=/home/siyuan/jobs/suzukii_WGS/rename_samples/cloned_input/sample_size.txt
conda_path=/home/siyuan/biosoft/miniconda3
conda_env=wgs_suzukii

bash /home/siyuan/jobs/suzukii_WGS/genetic_diff/ML_tree/results/mincount_10_no_se_auto_X_over10/TreeMix.sh \
-w $wd \
-s $dir_snp \
-c $chr_assignment \
-a $sample_size \
-l $window_len \
-p $script \
-d $conda_path \
-e $conda_env \
> run_winlen_$window_len\_mincov_$mincov\_mincount_$min_count\.log 2>&1
}

run_paratest 12 1 500
run_paratest 12 10 500
# run_paratest 12 5 500
# run_paratest 12 1 50
# run_paratest 12 1 75
# run_paratest 12 1 125