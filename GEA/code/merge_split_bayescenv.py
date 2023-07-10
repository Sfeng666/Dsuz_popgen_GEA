from optparse import OptionParser
import os
from statistics import mean
import time

### help and usage ### 
usage = "usage: %prog [options] args"
description = '''Function: merges EAA results from split BayeScEnv test jobs'''
version = '%prog 7.17.2022'
parser = OptionParser(usage=usage,version=version, description = description)                                 
(options,args) = parser.parse_args()

### variables ###
## paths ##
wd='/home/sfeng77/jobs/EAA/bayescenv/chtc/split_jobs/test_para_p/merge_jobs/results' # working directory
in_fst_dir='/home/sfeng77/jobs/EAA/bayescenv/chtc/split_jobs/test_para_p/results/{0}_{1}/out/{0}_{2}_{3}_fst.txt' # directory of output files of all split jobs
in_pos='/home/sfeng77/jobs/EAA/data/genotype_data/split/{}_pos_snp_gt_maf_5.txt' # position index of sites on autosomal or X-linked contigs

out_fst_signif='{0}/{1}_{2}/{1}_{3}_{4}_fst_signif_qval0.05.txt' # associated sites significant at qvalue < 0.05 in the same format
out_fst_signif_gsea='{0}/{1}_{2}/fst_signif_qval0.05_{1}_{3}.txt' # associated sites significant at qvalue < 0.05 in the GSEA input format
out_fst_signif_contig='{0}/{1}_{2}/contigs_fst_signif_qval0.05_{1}.txt' # contigs that have at least one significant associated site at qvalue < 0.05

out_log='{}/log.txt'.format(wd)
localtime = time.asctime( time.localtime(time.time()) )
os.system("echo 'Job starts at {0}' > {1}".format(localtime, out_log))

## list & dic ##
env_vars = ['Altitude', \
    'Annual_mean_wind_speed', \
    'Annual_Precipitation', 'Precipitation_Seasonality', \
    'Mean_Diurnal_Range', 'Mean_Temperature_of_Coldest_Quarter', 'Mean_Temperature_of_Warmest_Quarter', \
    'Ratio_Built_to_Vegetation']
# pr_prefs = [0.8, 0.5]
pr_prefs = [0.9]
chrs = ['auto', 'X']
genomic_pos = {chr: {} for chr in chrs}
contig2sites = {chr: {} for chr in chrs}
split_sites = {env_var: {} for env_var in env_vars}
split_sites_qval = {env_var: {} for env_var in env_vars}
split_run_outlier = {env_var: {} for env_var in env_vars} # check variation in # of outliers in each run
split_run_total = {env_var: {} for env_var in env_vars} # check variation in # of outliers in each run

with open(out_log, 'w') as f_log:
    header = '\t'.join(['environment_variable', '#_of_split', '#_of_outlier', '#_of_sampled_sites']) + '\n'
    f_log.write(header)

    ### 1. create output directions for each environmental variable
    for env_var in env_vars:
        for pr_pref in pr_prefs:
            os.system('mkdir -p {0}/{1}_{2}'.format(wd, env_var, pr_pref))

    ### 2. read the position index file for autosomal and X-linked contigs into a dictionary
    for chr in chrs:
        with open(in_pos.format(chr), 'r') as f_pos: # read the position index file for autosomal and X-linked contigs
            for line in f_pos:
                line = line.strip().split('\t')
                split_num = line[0]
                site_num = line[1]
                contig = line[2].split(':')[0]
                contig_pos = line[2].split(':')[1]
                genomic_pos[chr].setdefault(split_num)
                contig2sites[chr].setdefault(contig, {})
                contig2sites[chr][contig][contig_pos] = '_'.join([chr, split_num, site_num])

    for env_var in env_vars:
        for pr_pref in pr_prefs:
            fst_signif_contig=set()
            ### 2. read the output fst file into a dictionary
            for chr in chrs:
                for split_num in genomic_pos[chr]:
                    split_id = '_'.join([chr, split_num])
                    split_run_outlier[env_var].setdefault(split_id, 0) # check variation in # of outliers in each run
                    split_run_total[env_var].setdefault(split_id, 0) # check variation in # of outliers in each run
                    with open(in_fst_dir.format(env_var, pr_pref, chr, split_num), 'r') as f:
                        for line in f:
                            line = line.strip().split()
                            if line[0] == 'PEP_g':
                                header = '\t'.join(['genomic_pos', 'qval_g', 'PEP_g', 'g', 'PEP_alpha', 'alpha', 'fst']) + '\n'
                            else:
                                site_num = line[0]
                                pep_g = line[1]
                                g = line[3]
                                pep_a = line[4]
                                a = line[6]
                                fst = line[7]
                                split_sites[env_var]['_'.join([chr, split_num, site_num])] = [pep_g, g, pep_a, a, fst]
                                split_run_total[env_var][split_id] += 1 # check variation in # of outliers in each run

            ### 3. rank pep_g across all sites of the whole genome and recalculate qvalue for each site
            split_sites_sorted = sorted(split_sites[env_var], key=lambda x: float(split_sites[env_var][x][0]))
            split_sites_sorted_pep = list(float(split_sites[env_var][x][0]) for x in split_sites_sorted) # this step significantly save time for extracting pep value
            for i in range(0, len(split_sites_sorted)):
                qval = sum(split_sites_sorted_pep[:i + 1])/(i + 1)
                split_sites_qval[env_var][split_sites_sorted[i]] = [str(qval)] + split_sites[env_var][split_sites_sorted[i]]
                split_run_outlier[env_var]['_'.join(split_sites_sorted[i].split('_')[:2])] += 1 # check variation in # of outliers in each run
                if qval > 0.05: # only interested in sites with qval <= 0.05
                    break

            ### 4. output lines to each contig
            for chr in chrs:
                for contig in contig2sites[chr]:
                    with open(out_fst_signif.format(wd, env_var, pr_pref, chr, contig), 'w') as f, open(out_fst_signif_gsea.format(wd, env_var, pr_pref, contig), 'w') as f_gsea:
                        f.write(header)
                        f_gsea.write(header)
                        contig_pos_sorted = sorted(contig2sites[chr][contig], key=lambda x: int(x)) # sort the output file by genomic position of sites
                        for contig_pos in contig_pos_sorted:
                            if contig2sites[chr][contig][contig_pos] in split_sites_qval[env_var]:
                                line = [contig_pos] + split_sites_qval[env_var][contig2sites[chr][contig][contig_pos]]
                                f.write('\t'.join(line) + '\n')
                                f_gsea.write('\t'.join(line) + '\n')
                                fst_signif_contig.add(contig)
            
            ### 5. output contig ID that have significant sites
            with open(out_fst_signif_contig.format(wd, env_var, pr_pref), 'w') as f:
                for contig in fst_signif_contig:
                    f.write(contig + '\n')

            ### 6. output the number of outlier sites for each split run
            for split_id in split_run_outlier[env_var]:
                line = '\t'.join([env_var, split_id, str(split_run_outlier[env_var][split_id]), str(split_run_total[env_var][split_id])]) + '\n'
                f_log.write(line)

localtime = time.asctime(time.localtime(time.time()))
os.system("echo 'Job ends at {0}' >> {1}".format(localtime, out_log))