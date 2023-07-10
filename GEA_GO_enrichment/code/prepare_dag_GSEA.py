# -*- coding: utf-8 -*-
################################################################################
##  *
##  *  Function: prepare dag file to run GSEA on CHTC for top outlier sites]
##  *  Writer:Siyuan Feng
##  *  Mail: siyuanfeng.bioinfo@gmail.com
##  *  Created: 04.01.2022
##  *  Last modified: 10.03.2022
################################################################################

import os

def write_dag(env_var):
    dir_chtc = '/home/sfeng77/jobs'
    wd = '{0}/EAA/downstream/GSEA/random_SNP/scripts/output_gene_list/top_outliers/iterate/GSEA_top_500_genes/more_permutation_reps/results/{1}'.format(dir_chtc, env_var)
    dag_input = '{0}/split_jobs_{1}.dag'.format(wd, env_var)
    outdir_home = '{}/out'.format(wd)
    dir_para_scripts = '{0}/EAA/downstream/GSEA/random_SNP/scripts/output_gene_list/top_outliers/iterate/GSEA_top_500_genes/more_permutation_reps/scripts'.format(dir_chtc)
    dir_scripts = '{0}/EAA/downstream/GSEA/random_SNP/scripts'.format(dir_chtc)
    script = '{0}/GO_SNP_parallel_genomic_permutation_modifiedsuz_{1}.pl'.format(dir_scripts, env_var)
    in_outlier_prefix = '{0}/EAA/downstream/GSEA/random_SNP/results/output_gene_list/top_outliers/iterate/snp_subset_top500/{1}/top_outlier_snps_for_500_genes_{1}_'.format(dir_chtc, env_var)
    in_outlier_transfer = '{0}/EAA/downstream/GSEA/random_SNP/results/output_gene_list/top_outliers/iterate/snp_subset_top500/{1}/'.format(dir_chtc, env_var)
    in_anno_transfer = '{0}/EAA/downstream/GSEA/genomic_annot/'.format(dir_chtc)
    in_GO_gene = '{0}/EAA/downstream/GSEA/GO_annot/dmel/GO_gene_cats_parents20210708_sorted.txt'.format(dir_chtc)
    in_GO_desc = '{0}/EAA/downstream/GSEA/GO_annot/dmel/GO_desc12_condensed_parents.txt'.format(dir_chtc)
    in_randomSNP = '{0}/EAA/downstream/GSEA/random_SNP/snp_set_100M.txt'.format(dir_chtc)

    os.system('mkdir -p {}'.format(outdir_home))
    with open(dag_input, 'w') as f: 
        # set configuration variables for this dag, to remove limits on job submissions
        # note!! the assumption to allow unlimited job submissions, is there should not be too many job submissions at the same time.
        f.write('CONFIG {0}/unlimited.config\n'.format(dir_scripts))

        # 3.1. write paralleling job nodes for: running bayescenv
        request_cpu = 1
        request_disk = 8 * 1024**2
        request_memory = 1 * 1024
        
        for rep in range(1, 1001):    
            job = 'JOB {0}_{1} {2}/run_GSEA_chtc.sub\n'.format(env_var, rep, dir_para_scripts)
            vars = '''VARS {0}_{1} outdir_home="{2}" dir_scripts="{3}" script="{4}" dir_para_scripts="{14}" \
            in_outlier_prefix="{5}" in_outlier_transfer="{6}" in_anno_transfer="{7}" in_GO_gene="{8}" in_GO_desc ="{9}" in_randomSNP="{13}"\
            request_disk="{10}" request_memory="{11}" request_cpus="{12}" \
            env_var="{0}" rep="{1}" job="$(JOB)"\n'''.format(   
            env_var, rep, outdir_home, dir_scripts, script,
            in_outlier_prefix, in_outlier_transfer, in_anno_transfer, in_GO_gene, in_GO_desc,
            request_disk, request_memory, request_cpu, in_randomSNP, dir_para_scripts)
            f.writelines([job, vars])

###  write dag input for each batch of jobs with an environmental variable
env_vars = ['Altitude', \
    'Annual_mean_wind_speed', \
    'Annual_Precipitation', 'Precipitation_Seasonality', \
    'Mean_Diurnal_Range', 'Mean_Temperature_of_Coldest_Quarter', 'Mean_Temperature_of_Warmest_Quarter', \
    'Ratio_Built_to_Vegetation', 'Ratio_Crop_to_Forest']

for env_var in env_vars:
        write_dag(env_var)