
#!/bin/bash

## 1. receive arguments from submit_jobs.sh
env_var=$1
# BatchNum=0
script=/home/siyuan/jobs/suzukii_WGS/EAA/downstream/GSEA/random_SNP/scripts/output_gene_list/all_trimmed_outliers/scripts/GO_SNP_parallel_genomic_permutation_modifiedsuz_$env_var.pl
in_outlier_prefix=/home/siyuan/jobs/suzukii_WGS/EAA/bayescenv/chtc/split_jobs/test_para_p/merge_jobs/results/$env_var\_0.9/fst_signif_qval0.05_$env_var\_
path_out=/home/siyuan/jobs/suzukii_WGS/EAA/downstream/GSEA/random_SNP/scripts/output_gene_list/all_trimmed_outliers/results

## 2. set paths
wd_abs=$(pwd)

PermutationReps=1000
# script=$(basename $script)
# in_outlier_prefix=$(basename $in_outlier_prefix)

# output=$(basename $in_outlier_prefix)\_GO$BatchNum.txt
log=$path_out/log_$env_var.txt

# sleep 10m # test input file transfer

## 3. run perl scro[t for GSEA
echo -e "### GSEA for $env_var outlier replicate $BatchNum starts at $(date +%F'  '%H:%M) ###\n" > $log

perl $script \
-i $in_outlier_prefix \
-r $PermutationReps \
-e $path_out/$env_var \
>> $log 2>&1

echo -e "### GSEA for $env_var outlier replicate $BatchNum ends at $(date +%F'  '%H:%M) ###\n" >> $log

# ## 6. Handling output
# # For log files (could be redirected stdout/stderr) that you want to check under home directory, 
# # do not remove them at working directory, so that HTCondor will transfer them back to your home directory.
# shopt -s extglob
# rm -rf !($output|$log|_condor_stdout|_condor_stderr)
# shopt -u extglob