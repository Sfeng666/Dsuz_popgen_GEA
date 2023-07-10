#!/bin/bash

function run_paratest(){
mincov=$1
min_count=$2
wd=/home/siyuan/jobs/suzukii_WGS/reannotate/syn_diversity/mincov_$mincov\_mincount_$min_count
dir_snp=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/call_snp_maxcov_99/mincov_$mincov\_mincount_$min_count
dir_mpileup=/home/siyuan/jobs/suzukii_WGS/06mpileup_rsync
script=/home/siyuan/jobs/suzukii_WGS/calc_diversity_same_sites_ts/updated_poolsnp/filter_mpileup_by_cov.py
extract_script=/home/siyuan/jobs/suzukii_WGS/reannotate/syn_diversity/extract_4_syn_sites.py
bedtools=/opt/miniconda3/envs/play_with_genome/bin/bedtools
in_anno_seq=/home/siyuan/jobs/suzukii_WGS/reannotate/for_reannotate.fa
sample_size=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/calc_diversity_maxcov_99/test_normalize/sample_size.txt
maxcov=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/calc_maxcov_99/max_cov.txt
chr_assignment=/home/siyuan/jobs/suzukii_WGS/calc_maxcov/amb/assignment_cor_amb.txt
conda_path=/opt/miniconda3
conda_env=WGS_analysis

bash calc_diversity_synonymous.sh \
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
> run_mincov_$mincov\_mincount_$min_count\.log 2>&1
}

# run_paratest 12 10
run_paratest 12 1
# run_paratest 12 5
# run_paratest 12 8
# run_paratest 12 3
# run_paratest 12 20
# run_paratest 12 50
# run_paratest 20 10