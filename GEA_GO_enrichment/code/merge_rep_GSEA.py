# -*- coding: utf-8 -*-
################################################################################
##  *
##  *  Function: merge output p values fron replicates of GSEA]
##  *  Author: Siyuan Feng
##  *  Mail: siyuanfeng.bioinfo@gmail.com
##  *  Created: 04.19.2022
##  *  Last modified: 10.03.2022
################################################################################

import os

def merge_rep(env_var):
    dir_chtc = '/home/sfeng77/jobs'
    wd = '{0}/EAA/downstream/GSEA/random_SNP/scripts/output_gene_list/top_outliers/iterate/GSEA_top_500_genes/more_permutation_reps/results/{1}'.format(dir_chtc, env_var)
    outdir_home = '{}/out'.format(wd)
    out_rep = '{0}/top_outlier_snps_for_500_genes_{1}_GO{2}.txt'
    out_merge = '{0}/merged_GSEA_{1}.txt'.format(wd, env_var)
    alpha = 0.05
    out_merge_signif = '{0}/merged_GSEA_{1}_signif{2}.txt'.format(wd, env_var, alpha)
    go_info = {}

    for rep in range(1, 1001):
        with open(out_rep.format(outdir_home, env_var, rep), 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                go_id = line[0]
                p_val = float(line[5])
                if go_id in go_info:
                    go_info[go_id][5] += p_val
                else:
                    go_info[go_id] = line
                    go_info[go_id][5] = p_val

    for go in go_info:
        go_info[go][5] = go_info[go][5]/1000    # p-value calculated as the average among all replicate runs
    sorted_go = sorted(go_info, key=lambda go: go_info[go][5]) # sort GO terms by p-value

    with open(out_merge, 'w') as f, open(out_merge_signif, 'w') as f_signif:
        for go in sorted_go:
            row = '\t'.join(list(str(x) for x in go_info[go])) + '\n'
            f.write(row)
            if go_info[go][5] <= alpha:
                f_signif.write(row)

###  merge output p values fron replicates of GSEA for each environmental variable
env_vars = ['Altitude', \
    'Annual_mean_wind_speed', \
    'Annual_Precipitation', 'Precipitation_Seasonality', \
    'Mean_Diurnal_Range', 'Mean_Temperature_of_Coldest_Quarter', 'Mean_Temperature_of_Warmest_Quarter', \
    'Ratio_Built_to_Vegetation', 'Ratio_Crop_to_Forest']

for env_var in env_vars:
    merge_rep(env_var)