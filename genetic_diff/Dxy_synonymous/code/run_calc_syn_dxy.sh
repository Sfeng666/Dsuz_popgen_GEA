#!/bin/bash

function run_paratest(){
mincov=$1
min_count=$2
wd=/home/siyuan/jobs/suzukii_WGS/genetic_diff/Dxy/syn_Dxy/results
dir_snp=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/call_snp_maxcov_99/mincov_$mincov\_mincount_$min_count
dir_syn_snp=/home/siyuan/jobs/suzukii_WGS/genetic_diff/Fst/syn_Fst/results/syn_SNP_data
dir_mpileup=/home/siyuan/jobs/suzukii_WGS/06mpileup_rsync
script=/raid10/siyuan/marula_storage/jobs/suzukii_WGS/calc_diversity_same_sites_ts/updated_poolsnp/filter_mpileup_by_cov.py
calc_script=/home/siyuan/jobs/suzukii_WGS/genetic_diff/Dxy/calc_dxy_contig.py
maxcov=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/calc_maxcov_99/max_cov.txt
chr_assignment=/home/siyuan/jobs/suzukii_WGS/calc_maxcov/amb/assignment_cor_amb.txt
conda_path=/opt/miniconda3
conda_env=WGS_analysis

bash $wd/../scripts/calc_syn_dxy.sh \
-w $wd \
-s $dir_snp \
-y $dir_syn_snp \
-m $dir_mpileup \
-c $chr_assignment \
-v $mincov \
-x $maxcov \
-t $min_count \
-l $calc_script \
-p $script \
-d $conda_path \
-e $conda_env \
> $wd/../run_calc_syn_dxy.log 2>&1
}

# run_paratest 12 10
run_paratest 12 1
# run_paratest 12 5
# run_paratest 12 8
# run_paratest 12 3
# run_paratest 12 20
# run_paratest 12 50
# run_paratest 20 10