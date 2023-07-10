#!/bin/bash

function run_paratest(){
mincov=$1
min_count=$2
wd=/home/siyuan/jobs/suzukii_WGS/genetic_diff/Fst/syn_Fst/results/syn_SNP_data
dir_snp=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/call_snp_maxcov_99/mincov_$mincov\_mincount_$min_count
dir_mpileup=/home/siyuan/jobs/suzukii_WGS/06mpileup_rsync
script=/home/siyuan/jobs/suzukii_WGS/calc_diversity_same_sites_ts/updated_poolsnp/filter_mpileup_by_cov.py
extract_script=/raid10/siyuan/marula_storage/jobs/suzukii_WGS/reannotate/syn_diversity/extract_4_syn_sites.py
bedtools=/opt/miniconda3/envs/play_with_genome/bin/bedtools
in_anno_seq=/raid10/siyuan/marula_storage/jobs/suzukii_WGS/reannotate/for_reannotate.fa
sample_size=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/calc_diversity_maxcov_99/test_normalize/sample_size.txt
maxcov=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/calc_maxcov_99/max_cov.txt
chr_assignment=/home/siyuan/jobs/suzukii_WGS/calc_maxcov/amb/assignment_cor_amb.txt
conda_path=/opt/miniconda3
conda_env=WGS_analysis

bash /home/siyuan/jobs/suzukii_WGS/genetic_diff/Fst/syn_Fst/scripts/extract_4_syn_sites.sh \
-w $wd \
-s $dir_snp \
-n $sample_size \
-m $dir_mpileup \
-c $chr_assignment \
-v $mincov \
-x $maxcov \
-t $min_count \
-a $in_anno_seq \
-r $extract_script \
-p $script \
-b $bedtools \
-d $conda_path \
-e $conda_env \
> $wd/../run_mincov_$mincov\_mincount_$min_count\.log 2>&1
}

run_paratest 12 1