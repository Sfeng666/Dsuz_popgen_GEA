#!/bin/bash

## 1. receive arguments from run_bayescenv_chtc.sub
chr=$1
env_var=$2
num=$3
in_snp=$4
in_env=$5
pr_pref=$6

## 2. set paths
wd_abs=$(pwd)
bayscnenv=biosoft/bayescenv-1.1/bin/linux64/bayescenv

prefix=$env_var\_$chr\_$num
out_stat=$prefix\_fst.txt
out_acc=$prefix\_AccRte.txt
out_sel=$prefix.sel
out_verif=$prefix\_Verif.txt

log=log_$prefix.txt

## 3. run BayeScEnv
echo -e "### BayeScEnv on $chr chromosome, subsample $num with $env_var starts at $(date +%F'  '%H:%M) ###\n" >> $log

$bayscnenv $(basename $in_snp) \
-pr_jump 0.02 \
-pr_pref $pr_pref \
-env $(basename $in_env) \
-od $wd_abs \
-o $prefix \
>> $log 2>&1
echo -e "### BayeScEnv on $chr chromosome, subsample $num with $env_var ends at $(date +%F'  '%H:%M) ###\n" >> $log

## 6. Handling output
# For log files (could be redirected stdout/stderr) that you want to check under home directory, 
# do not remove them at working directory, so that HTCondor will transfer them back to your home directory.
shopt -s extglob
rm -rf !($out_stat|$out_acc|$out_sel|$out_verif|$log)
shopt -u extglob