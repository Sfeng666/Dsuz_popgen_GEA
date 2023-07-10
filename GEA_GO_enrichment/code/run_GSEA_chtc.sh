#!/bin/bash

## 1. receive arguments from run_bayescenv_chtc.sub
env_var=$1
BatchNum=$2
script=$3
in_outlier_prefix=$4

## 2. set paths
wd_abs=$(pwd)

# PermutationReps=1000
PermutationReps=100000 # try more permutation reps to get a more precise P value
script=$(basename $script)
in_outlier_prefix=$(basename $in_outlier_prefix)

output=$(basename $in_outlier_prefix)\_GO$BatchNum.txt
log=log_$env_var\_$BatchNum.txt

# sleep 10m # test input file transfer

## 3. run perl script for GSEA
echo -e "### GSEA for $env_var outlier replicate $BatchNum starts at $(date +%F'  '%H:%M) ###\n" >> $log

perl $script \
-i $in_outlier_prefix \
-b $BatchNum \
-r $PermutationReps \
>> $log 2>&1

echo -e "### GSEA for $env_var outlier replicate $BatchNum ends at $(date +%F'  '%H:%M) ###\n" >> $log
rm $log

# ## 6. Handling output
# # For log files (could be redirected stdout/stderr) that you want to check under home directory, 
# # do not remove them at working directory, so that HTCondor will transfer them back to your home directory.
# shopt -s extglob
# rm -rf !($output)
# shopt -u extglob