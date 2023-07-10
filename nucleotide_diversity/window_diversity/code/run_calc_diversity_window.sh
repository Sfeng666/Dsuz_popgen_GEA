#!/bin/bash

function run_paratest(){
mincov=$1
min_count=$2
win_len=$3
wd=/home/siyuan/jobs/suzukii_WGS/genetic_diff/window_diversity/results/winlen_$win_len\_mincov_$mincov\_mincount_$min_count
sd=/home/siyuan/jobs/suzukii_WGS/genetic_diff/window_diversity/scripts
script=$sd/win_div_by_filtered_sites.py
calc_script=$sd/calc_diversity_window.py
dir_snp=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/call_snp_maxcov_99/mincov_$mincov\_mincount_$min_count
dir_mpileup=/home/siyuan/jobs/suzukii_WGS/06mpileup_rsync
chr_assignment=/home/siyuan/jobs/suzukii_WGS/calc_maxcov/amb/assignment_cor_amb.txt
sample_size=/home/siyuan/jobs/suzukii_WGS/rename_samples/cloned_input/sample_size.txt
in_assembly_rep=/home/siyuan/reference/WT3-2.0/refseq/annotation/GCF_013340165.1_LBDM_Dsuz_2.1.pri_assembly_report.txt
maxcov=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/calc_maxcov_99/max_cov.txt
conda_path=/opt/miniconda3
conda_env=WGS_analysis

bash $sd/calc_diversity_window.sh \
-w $wd \
-s $dir_snp \
-m $dir_mpileup \
-c $chr_assignment \
-a $sample_size \
-r $in_assembly_rep \
-f $win_len \
-v $mincov \
-x $maxcov \
-t $min_count \
-l $calc_script \
-p $script \
-d $conda_path \
-e $conda_env \
> run_winlen_$win_len\_mincov_$mincov\_mincount_$min_count\.log 2>&1
}

run_paratest 12 1 125000
# run_paratest 12 1 150
# run_paratest 12 1 200
# run_paratest 12 1 50
# run_paratest 12 1 75
# run_paratest 12 1 125