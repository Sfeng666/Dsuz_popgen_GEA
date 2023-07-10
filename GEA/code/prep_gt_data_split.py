from optparse import OptionParser
import math
import sympy

### help and usage ### 
usage = "usage: %prog [options] args"
description = '''Function: prepare the split input genotype data for BayeScEnv'''
version = '%prog 2.28.2022'
parser = OptionParser(usage=usage,version=version, description = description)
parser.add_option("--wd",
                    action="store",
                    dest = "wd",
                    help = "working directory",
                    metavar = "PATH")
parser.add_option("--in_snp_dir",
                    action="store",
                    dest = "in_snp_dir",
                    help = "the directory of filtered SNPs from all contigs in vcf format (.vcf)",
                    metavar = "PATH")
parser.add_option("--in_sample_size",
                    action="store",
                    dest = "in_sample_size",
                    help = "input path of allelic (auto/X) sample sizes of all samples analyzed (.txt)",
                    metavar = "PATH")
parser.add_option("--in_chr_assign",
                    action="store",
                    dest = "in_chr_assign",
                    help = "input path of assignments of contigs to chromosomes (.txt)",
                    metavar = "PATH")
parser.add_option("--in_sample_order",
                    action="store",
                    dest = "in_sample_order",
                    help = "input path of a file that contains the order of samples as environmental variables (.txt)",
                    metavar = "PATH")
parser.add_option("--out_log",
                    action="store",
                    dest = "out_log",
                    help = "output path of the log file (.log)",
                    metavar = "PATH")                                        
(options,args) = parser.parse_args()

### introduced variables ###
wd = options.wd
in_snp_dir = options.in_snp_dir
in_sample_size = options.in_sample_size
in_chr_assign = options.in_chr_assign
in_sample_order = options.in_sample_order
out_log = options.out_log

### local variables ###
subsample_size = 10000

### def functions ###
## function to calculate effective sample size as the expectation of unique chromosomal draws
def expectation_sympy(nr, nc):
    return float(sum(j*math.factorial(nc)*sympy.functions.combinatorial.numbers.stirling(nr, j)/(math.factorial(nc - j)*nc**nr) for j in range(1, nc + 1)))

### 1. input chromosomal assignment information of contigs
chrs = ['X', 'auto']
chr_assign = {x:{} for x in chrs}
with open(in_chr_assign, 'r') as f:
    for line in f:
        if not line.startswith('contig'):
            line = line.strip().split('\t')
            contig = line[0]
            mapping = line[1]
            cor_assign = line[2]
            assign = cor_assign.split('-')[0]
            if mapping != 'na':
                if not mapping in chr_assign[assign]:
                    chr_assign[assign][mapping] = [contig]
                else:
                    chr_assign[assign][mapping].append(contig)
            else:
                if cor_assign == 'X-linked':
                    if not 'X' in chr_assign[assign]:
                        chr_assign[assign]['X'] = [contig]
                    else:
                        chr_assign[assign]['X'].append(contig)
                elif cor_assign == 'auto-linked':
                    if not cor_assign in chr_assign[assign]:
                        chr_assign[assign][cor_assign] = [contig]
                    else:
                        chr_assign[assign][cor_assign].append(contig)

### 2. extract haplotype size for each sample
sample_size = {}
with open(in_sample_size, 'r') as f:
    i = 0
    for line in f:
        line = line.strip().split('\t')
        if i == 0:
            samples = line
        elif i == 1:
            sample_size['auto'] = {x: int(y) for x,y in zip(samples, line)}
        elif i == 2:
            sample_size['X'] = {x: int(y) for x,y in zip(samples, line)}
        i += 1

### 3. extract the order of samples in the input of environmental variables
sample_order = []
with open(in_sample_order, 'r') as f:
    for line in f:
        if not line.startswith('Population'):
            line = line.strip().split('\t')
            sample_order.append(line[0])

### 4. generate genotype file for BayeScEnv
expectation_dict = {}
for chr in chrs:
    dic_gt = {sample: {} for sample in sample_order}
    num_snp, num_snp_filter = 0, 0
    out_gt_pos = '{0}/{1}_pos_snp_gt_maf_5.txt'.format(wd, chr)
    for arm in chr_assign[chr]:
        for contig in chr_assign[chr][arm]:
            in_snp = '{0}/{1}_snp_clean.vcf'.format(in_snp_dir, contig)
            with open(in_snp, 'r') as f:
                for line in f:
                    if line.startswith('#CHROM'):
                        samples = line.strip().split(' ')[-1].split('\t')
                    elif not line.startswith('#'):
                        line = line.strip().split('\t')
                        num_snp += 1
                        pos = int(line[1])
                        ref = line[3].upper()
                        alt = line[4].upper()
                        alleles = [ref] + alt.split(',')
                        allele_freq = {x:0 for x in alleles}
                        ct_sp = {}

                        # extract allele frequency of the top 2 alleles across all samples
                        for i in range(len(line[9:])):
                            info = line[9:][i]
                            gt = list(int(x) for x in info.split(':')[0].split('/'))
                            ct = [int(info.split(':')[1])] + list(int(x) for x in info.split(':')[2].split(','))
                            dp = int(info.split(':')[3])

                            if len(set(gt)) == 1:
                                ct_sp[i] = {alleles[gt[0]]: max(ct)}
                            else:
                                ct_sp[i] = {alleles[x]: ct[gt.index(x)] for x in gt}

                            for ale in alleles:
                                if ale in ct_sp[i]:
                                    allele_freq[ale] += ct_sp[i][ale]/dp
                                else:
                                    ct_sp[i][ale] = 0

                        # sort alleles by the total allele frequency (so that each sample has the same weight of voting despite sequencing depth), and determine the minor allele
                        bialleles = sorted(allele_freq, key=lambda x: allele_freq[x], reverse=True)[:2]

                        # extract allele count (of top two alleles) of each sample, corrected by effective sample size
                        cts_maj = []
                        cts_min = []
                        sum_min_ct = 0
                        sum_bi_ct = 0
                        for i in range(len(line[9:])):
                            dp_bi = sum(list(ct_sp[i][ale] for ale in bialleles if ale in ct_sp[i]))
                            sum_bi_ct += dp_bi
                            if dp_bi != 0:
                                af_major = ct_sp[i][bialleles[0]]/dp_bi
                                af_minor = ct_sp[i][bialleles[1]]/dp_bi
                                sum_min_ct += ct_sp[i][bialleles[1]]
                                
                                # only calculate effective sample size for combinations of depth and pool size that haven't been calculated
                                key_sz = ','.join([str(dp_bi), str(sample_size[chr][samples[i]])])
                                if key_sz in expectation_dict:
                                    efs = expectation_dict[key_sz]
                                else:
                                    efs = expectation_sympy(dp_bi, sample_size[chr][samples[i]])
                                    expectation_dict[key_sz] = efs
                                ct_major = af_major*efs
                                ct_minor = af_minor*efs

                            else:
                                ct_major, ct_minor = 0, 0
                            cts_maj.append(ct_major)
                            cts_min.append(ct_minor)
                        
                        # filter SNPs by across-sample MAF >= 5%
                        if sum_min_ct/sum_bi_ct >= 0.05:
                            num_snp_filter += 1
                            for sample in dic_gt:
                                dic_gt[sample][len(dic_gt[sample]) + 1] = [':'.join([contig, str(pos)]), \
                                    str(round(cts_maj[samples.index(sample)]) + round(cts_min[samples.index(sample)])), \
                                    '2', str(round(cts_maj[samples.index(sample)])), \
                                        str(round(cts_min[samples.index(sample)]))]
    ## subsample SNPs by a fixed interval
    with open(out_gt_pos, 'w') as f_pos, open(out_log, 'a') as f_log:
        # output the filtering statictics
        f_log.write('{0}: Number of SNPs = {1}\n'.format(chr, num_snp))
        f_log.write('{0}: Number of SNPs after MAF filtering = {1}\n'.format(chr, num_snp_filter))
        f_log.write('{0}: Proportion of SNPs filtered out = {1:.2%}\n'.format(chr, (num_snp - num_snp_filter)/num_snp))

        # determine the interval length of subsampling by sample size. interval length = numebr of samples
        interval = math.ceil(len(dic_gt[sample_order[0]])/subsample_size)
        subsamplesize_div = math.ceil(len(dic_gt[sample_order[0]])/interval)
        f_log.write('{0}: Number of subsamples = {1}\n'.format(chr, interval))

        # output subsampled SNPs
        for sub in range(1, interval + 1):
            out_gt = '{0}/{1}_{2}_snp_gt_maf_5.txt'.format(wd, chr, sub)
            with open(out_gt, 'w') as f:
                if sub + (subsamplesize_div - 1)*interval <= len(dic_gt[sample_order[0]]):
                    subsamplesize_actual = subsamplesize_div
                else:
                    subsamplesize_actual = subsamplesize_div - 1
                f.write('[loci]={}\n\n'.format(subsamplesize_actual))
                f.write('[populations]={}\n\n'.format(len(sample_order)))
                for pop in sample_order:
                    f.write('[pop]={}\n'.format(sample_order.index(pop) + 1))
                    for i in range(subsamplesize_actual):
                        f.write('\t'.join([str(i + 1)] + dic_gt[pop][sub + i*interval][1:]) + '\n')
                        if pop == sample_order[0]: # avoid repeated output to the position file
                            f_pos.write(('\t'.join([str(sub), str(i + 1), dic_gt[pop][sub + i*interval][0]]) + '\n'))
                    f.write('\n')
