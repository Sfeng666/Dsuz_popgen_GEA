# script to replace values of the @chrs, @StartSites, @StopSites perl variables with suzukii contig information

import os
import re

env_vars = ['Altitude', \
    'Annual_mean_wind_speed', \
    'Annual_Precipitation', 'Precipitation_Seasonality', \
    'Mean_Diurnal_Range', 'Mean_Temperature_of_Coldest_Quarter', 'Mean_Temperature_of_Warmest_Quarter', \
    'Ratio_Built_to_Vegetation', 'Ratio_Crop_to_Forest']
# contig_signif='/home/siyuan/jobs/suzukii_WGS/EAA/bayescenv/chtc/merge_jobs_corrected/results/{0}/contigs_fst_signif_qval0.05_{0}.txt'
dir_GEA = '/home/siyuan/jobs/suzukii_WGS/EAA/bayescenv/chtc/merge_jobs_corrected/results'
in_gff='/home/siyuan/reference/WT3-2.0/refseq/annotation/GCF_013340165.1_LBDM_Dsuz_2.1.pri_genomic.gff'
in_script='/home/siyuan/jobs/suzukii_WGS/EAA/downstream/GSEA/random_SNP/scripts/GO_SNP_parallel_genomic_permutation_modifiedsuz.pl'
out_script='/home/siyuan/jobs/suzukii_WGS/EAA/downstream/GSEA/random_SNP/scripts/GO_SNP_parallel_genomic_permutation_modifiedsuz_{0}.pl'

# 1. get the name of contigs that contain analyzed SNP
for env_var in env_vars:
    contigs = {}
    cmd = "echo {0}/{1}/fst_signif_qval0.05_{1}_*".format(dir_GEA, env_var)
    file_contigs = os.popen(cmd).read().strip().split()
    for file_contig in file_contigs:
        contig = "'" + re.search('fst_signif_qval0.05_{0}_(.*).txt'.format(env_var), file_contig).group(1) + "'"
        contigs[contig] = []         

    with open(in_gff, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                contig = "'" + line[0] + "'"
                type = line[2]
                end = line[4]
                if type == 'region' and contig in contigs:
                    contigs[contig] += ['1', end]

    os.system('''sed "s/'Chr2L','Chr2R','Chr3L','Chr3R','ChrX'/{0}/g" {3} > {4}
    sed -i  "s/1,1,1,1,1/{1}/g" {4}
    sed -i  "s/23011543,21146707,24543556,27905052,22422826/{2}/g" {4}'''.format(','.join(contigs), 
    ','.join(list(contigs[contig][0] for contig in contigs)), 
    ','.join(list(contigs[contig][1] for contig in contigs)),
    in_script, out_script.format(env_var)))

