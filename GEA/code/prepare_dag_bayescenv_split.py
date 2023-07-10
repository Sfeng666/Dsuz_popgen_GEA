# -*- coding: utf-8 -*-
################################################################################
##  *
##  *  Function: prepare dag file to run split BayeScEnv on CHTC]
##  *  Writer:Siyuan Feng
##  *  Mail: siyuanfeng.bioinfo@gmail.com
##  *  Version: 03.01.2021
##  *
################################################################################

import os
import re

def write_dag(env_var, cpu, longjob, pr_pref):
    wd = '/home/sfeng77/jobs/EAA/bayescenv/chtc/split_jobs/test_para_p/results/{0}_{1}'.format(env_var, pr_pref)
    dag_input = '{0}/split_jobs_{1}_{2}.dag'.format(wd, env_var, pr_pref)
    outdir_home = '{}/out'.format(wd)
    dir_scripts = '/home/sfeng77/jobs/EAA/bayescenv/chtc/split_jobs/test_para_p/scripts'
    software_pack = '/home/sfeng77/biosoft/biosoft'
    in_snp_pre = '/home/sfeng77/jobs/EAA/data/genotype_data/split'
    in_env_pre = '/home/sfeng77/jobs/EAA/data/env_data/input_BayeScEnv'
    in_snp_suf ='_snp_gt_maf_5.txt'
    in_env_suf = '.txt'
    
    os.system('mkdir -p {}'.format(outdir_home))
    with open(dag_input, 'w') as f: 
        # set configuration variables for this dag, to remove limits on job submissions
        # note!! the assumption to allow unlimited job submissions, is there should not be too many job submissions at the same time.
        f.write('CONFIG {0}/unlimited.config\n'.format(dir_scripts))

        # 3.1. write paralleling job nodes for: running bayescenv
        request_disk = 1 * 1024**2
        request_memory = 1.5 * 1024
        
        for chr in sp_num:
            for num in range(1, sp_num[chr] + 1):    
                in_snp = '{0}/{1}_{2}{3}'.format(in_snp_pre, chr, num, in_snp_suf)
                in_env = '{0}/{1}{2}'.format(in_env_pre, env_var, in_env_suf)
                job = 'JOB {0}_{1}_{2} {3}/run_bayescenv_chtc.sub\n'.format(env_var, chr, num, dir_scripts)
                vars = '''VARS {0}_{1}_{2} chr="{1}" env_var="{0}" in_snp="{3}" in_env="{4}" \
                outdir_home="{5}" request_disk="{6}" request_memory="{7}" request_cpus="{8}" \
                software_pack="{9}" dir_scripts="{10}" longjob="{11}" pr_pref="{12}" num="{2}" \
                job="$(JOB)"\n'''.format(
                env_var, chr, num, in_snp, in_env, outdir_home, request_disk, request_memory, 
                cpu, software_pack, dir_scripts, longjob, pr_pref)
                f.writelines([job, vars])

### 1. Define the list of environmental variables and the list of CPUs to run each environmenta batch  ###
env_vars = ['Altitude', \
    'Annual_mean_wind_speed', \
    'Annual_Precipitation', 'Precipitation_Seasonality', \
    'Mean_Diurnal_Range', 'Mean_Temperature_of_Coldest_Quarter', 'Mean_Temperature_of_Warmest_Quarter', \
    'Ratio_Built_to_Vegetation']
cpus = [1]*len(env_vars)
longjobs = ['true']*len(env_vars)
pr_prefs = [0.9, 0.8]

### 2. extract the number of subsamples for X and autosomes
sp_num_log = '/home/sfeng77/jobs/EAA/data/genotype_data/split/log.txt'
sp_num = {}
with open(sp_num_log, 'r') as f:
    for line in f:
        if re.search('auto: Number of subsamples = ', line):
            sp_num['auto'] = int(re.search('auto: Number of subsamples =(.*$)', line).group(1))
        elif re.search('X: Number of subsamples = ', line):
            sp_num['X'] = int(re.search('X: Number of subsamples =(.*$)', line).group(1))

### 3. write dag input for each batch of jobs with an environmental variable
for env_var, cpu, longjob in zip(env_vars, cpus, longjobs):
    for pr_pref in pr_prefs:
        write_dag(env_var, cpu, longjob, pr_pref)