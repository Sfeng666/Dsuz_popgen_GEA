# script to generate tables of full information about significantly associated loci and genes

import os
import re
import xlsxwriter

# file paths and directories
dir_GEA = '/home/siyuan/jobs/suzukii_WGS/EAA/bayescenv/chtc/split_jobs/test_para_p/merge_jobs/results'
dir_snp_gene='/home/siyuan/jobs/suzukii_WGS/EAA/downstream/GSEA/random_SNP/scripts/output_gene_list/all_trimmed_outliers/results'
script='/home/siyuan/jobs/suzukii_WGS/EAA/downstream/GSEA/random_SNP/scripts/output_gene_list/top_outliers/iterate/run_GSEA_top10.sh'
in_ortho = '/home/siyuan/reference/WT3-2.0/refseq/ortholog/gene_ortholog_suzukii_melanogaster.txt'
in_chr_assign = '/home/siyuan/jobs/suzukii_WGS/calc_maxcov/amb/assignment_cor_amb.txt'

# initial variables
env_vars = ['Altitude', \
    'Annual_mean_wind_speed', \
    'Annual_Precipitation', 'Precipitation_Seasonality', \
    'Mean_Diurnal_Range', 'Mean_Temperature_of_Coldest_Quarter', 'Mean_Temperature_of_Warmest_Quarter', \
    'Ratio_Built_to_Vegetation', 'Ratio_Crop_to_Forest']
abbrs = ["Altitude", "Mean wind speed", "Annual precip.", "Precip. seasonality", 
                  "Mean diurnal range", "Mean temp. cold qtr.", "Mean temp. warm qtr.", 
                  "Ratio built to vegetation", "Ratio crop to forest"]  # ":" is not allowed in work sheet name, so replaced with 'to'
env_abbrs = {env : abbr for env, abbr in zip(env_vars, abbrs)}

# 1. read the orthology table
ortho_fb = {}   # dictionary of Dmel - Dsuz gene orthology, with FBgn as keys
with open(in_ortho, 'r') as f:
    for line in f:
            line = line.strip().split('\t')
            gene = line[0]
            symb = line[2]
            fb = line[4].split(':')[1]
            des = line[5]
            ortho_fb[fb] = [gene, symb, des]

# 2. input chromosomal assignment information of contigs, in the form of {contig: [assigned_chromosome(auto/X), mapped_chromosome(2R/2L/3R/3L/X/na)]}
contig_assign = {}
with open(in_chr_assign, 'r') as f:
    for line in f:
        if not line.startswith('contig'):
            line = line.strip().split('\t')
            contig = line[0]
            arm = line[1]
            chr = line[2].split('-')[0]
            contig_assign[contig] = [chr, arm]

workbook = xlsxwriter.Workbook('{0}/sup_table_all_envs_snp_gene_stats.xlsx'.format(dir_snp_gene)) # create an excel table for all environnmental variable
for env in env_vars:
    worksheet = workbook.add_worksheet(env_abbrs[env]) # create a worksheet for each environmental variable
    # 3. extract stats of outlier sites from all contigs
    snp_contig = {}
    dic_g = {}
    cmd = "echo {0}/{1}_0.9/fst_signif_qval0.05_{1}_*.txt".format(dir_GEA, env)
    file_contigs = os.popen(cmd).read().strip().split()
    for file_contig in file_contigs:
        contig = re.search('fst_signif_qval0.05_{0}_(.*).txt'.format(env), file_contig).group(1)
        with open(file_contig, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                if line[0] != 'genomic_pos':
                    pos = line[0]
                    qval_g = float(line[1])
                    g = float(line[3])
                    snp_contig[':'.join([contig, pos])] = qval_g
                    dic_g[':'.join([contig, pos])] = g
    
    # 4. extract snp-gene association for trimmed outliers
    file_snp_gene = '{0}/{1}_SNP_gene_table.txt'.format(dir_snp_gene, env)
    snp_gene = {}
    with open(file_snp_gene, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            contig  = line[0]
            pos = line[1]
            if len(line) != 2:  # check if there are genes associated with this site
                genes = list(ortho_fb[fb][1] for fb in line[2].split()) # transfer FBgns to gene symbol
            else:
                genes = []
            snp_gene[':'.join([contig, pos])] = genes

    # 5. sort trimmed outliers first by qvalue_g (increasing), then by g (decreasing)
    sorted_snps = sorted(snp_gene, key = lambda x: (snp_contig[x], -dic_g[x]))

    # 6. output a table containing six columns: chromosome/chromosomal arm, contig, SNP position, qvalue_g, g, associated genes.
    file_out = '{0}/{1}_SNP_gene_stats.txt'.format(dir_snp_gene, env)
    with open(file_out, 'w') as f:
        header = ['Chromosome', 'Contig', 'Position', 'Q-value of g parameter', 'G parameter', 'Associated genes']
        f.write('\t'.join(header) + '\n')
        row = 0
        col = 0
        for item in header:
            worksheet.write(row, col, item)
            col += 1
        row = 1
        for snp in sorted_snps:
            contig = snp.split(':')[0]
            pos = snp.split(':')[1]
            if contig_assign[contig][1] != 'na':    # output either chromosomal arm (if directly mapped) or chromosome (auto or X) for plotting
                chr = contig_assign[contig][1]
            else:
                chr = contig_assign[contig][0] .replace('auto', 'autosome')
            list_row = [chr, contig, int(pos), snp_contig[snp], dic_g[snp], ','.join(snp_gene[snp])]
            f.write('\t'.join(list(str(x) for x in list_row)) + '\n')
            col = 0
            for item in list_row:
                worksheet.write(row, col, item)
                col += 1
            row += 1
workbook.close()