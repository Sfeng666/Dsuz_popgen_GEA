from optparse import OptionParser
from math import factorial
import sympy

### help and usage ### 
usage = "usage: %prog [options] args"
description = '''Function: prepare the input table of SNP allele frequency for TreeMix'''
version = '%prog 09.29.2021'
parser = OptionParser(usage=usage,version=version, description = description)
parser.add_option("--in_snp",
                    action="store",
                    dest = "in_snp",
                    help = "input path of filtered SNPs from a specified chr in vcf format (.vcf)",
                    metavar = "PATH")
parser.add_option("--in_sample_size",
                    action="store",
                    dest = "in_sample_size",
                    help = "input path of allelic (auto/X) sample sizes of all samples analyzed (.txt)",
                    metavar = "PATH")
parser.add_option("--out_window_len",
                    action="store",
                    dest = "out_window_len",
                    help = "output path of the table of window length of all contigs (.txt)",
                    metavar = "PATH")
parser.add_option("--out_SNP_num",
                    action="store",
                    dest = "out_SNP_num",
                    help = "output path of the file recording number of SNPs within each contig (.txt)",
                    metavar = "PATH")
parser.add_option("--chr",
                    action="store",
                    dest = "chr",
                    help = "a string indicating whether the contig being processed is auto- or X-linked",
                    metavar = "STR")
parser.add_option("--temp",
                    action = "store",
                    dest = "temp",
                    help = "output directory of temporary intermediate files",
                    metavar = "PATH")
(options,args) = parser.parse_args()

### introduced variables ###
in_snp = options.in_snp
in_sample_size = options.in_sample_size
out_window_len = options.out_window_len
out_SNP_num = options.out_SNP_num
chr = options.chr
temp = options.temp

### def functions ###
## function to calculate effective sample size as the expectation of unique chromosomal draws
def expectation_sympy(nr, nc):
    return float(sum(j*factorial(nc)*sympy.functions.combinatorial.numbers.stirling(nr, j)/(factorial(nc - j)*nc**nr) for j in range(1, nc + 1)))

### 1. read the input of sample size (the total number of chromosomes in a pool)
sample_size = {'auto': [], 'X': []}
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
        
### 2. count the window length based on a set of number of SNPs included in each window
set_snpnum = range(1, 501) 
start, end, num_snp, win = {x: 0 for x in set_snpnum}, {x: 0 for x in set_snpnum}, {x: 0 for x in set_snpnum}, {x: 0 for x in set_snpnum}
windows = {x: {0: {'start': 1}} for x in set_snpnum}
num_snp_contig = 0

expectation_dict = {}
first = 1
with open(in_snp, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            num_snp = {x: num_snp[x] + 1 for x in set_snpnum}
            num_snp_contig += 1
            pos = int(line[1])
            ref = line[3].upper()
            alt = line[4].upper()
            alleles = [ref] + alt.split(',')
            allele_freq = {x:0 for x in alleles}
            ct_sp = {}

            if first == 1:
                windows = {x: {0: {'start': pos}} for x in set_snpnum}
                first = 0

            # specify the ending boundry of each window at the midpoint of the end SNP and the next starting SNP
            # specify the starting boundry of each window at the midpoint of the starting SNP and the previous ending SNP + 1
            for num_cutoff in set_snpnum:
                if end[num_cutoff] == 1 and start[num_cutoff] == 1:
                    sum_pos = windows[num_cutoff][win[num_cutoff]]['end'] + pos
                    if sum_pos & 1 == 0:
                        windows[num_cutoff][win[num_cutoff]]['end'] = int(sum_pos/2)
                    else:
                        windows[num_cutoff][win[num_cutoff]]['end'] = int((sum_pos - 1)/2)
                    win[num_cutoff] += 1
                    windows[num_cutoff].setdefault(win[num_cutoff], {'start': windows[num_cutoff][win[num_cutoff] - 1]['end'] + 1})
                    start[num_cutoff] = 0
                    end[num_cutoff] = 0

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

            # sort alleles by the total allele frequency, and determine the minor allele
            bialleles = sorted(allele_freq, key=lambda x: allele_freq[x], reverse=True)[:2]

            # extract allele count (of top two alleles) of each sample, corrected by effective sample size
            cts = []
            for i in range(len(line[9:])):
                dp_bi = sum(list(ct_sp[i][ale] for ale in bialleles if ale in ct_sp[i]))
                if dp_bi != 0:
                    af_major = ct_sp[i][bialleles[0]]/dp_bi
                    if bialleles[1] in ct_sp[i]:
                        af_minor = ct_sp[i][bialleles[1]]/dp_bi
                    else:
                        af_minor = 0
                    
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
                cts.append('{0},{1}'.format('{:.2f}'.format(ct_major), '{:.2f}'.format(ct_minor)))
            # f_out.write(' '.join([in_snp] + cts) + '\n')
            # print(' '.join([in_snp] + cts))
            print(' '.join(cts))

            # partition windows based on the number of SNPs contained in each window
            for num_cutoff in set_snpnum:
                if num_snp[num_cutoff] == num_cutoff:
                    windows[num_cutoff][win[num_cutoff]]['end'] = pos
                    windows[num_cutoff][win[num_cutoff]]['num_snp'] = num_snp[num_cutoff]
                    end[num_cutoff] = 1
                    start[num_cutoff] = 1

                    num_snp[num_cutoff] = 0

    for num_cutoff in set_snpnum:
        if num_snp[num_cutoff] < num_cutoff:
            if num_snp[num_cutoff] >= num_cutoff*0.95:
                windows[num_cutoff][win[num_cutoff]]['end'] = pos
                windows[num_cutoff][win[num_cutoff]]['num_snp'] = num_snp[num_cutoff]
            else:
                del windows[num_cutoff][win[num_cutoff]]

# with open(out_SNP_num, 'a') as f:
#     f.write(' '.join([in_snp, str(num_snp_contig)]) + '\n')

### 4. output the length of windows
with open(out_window_len, 'a') as f:
    for num_cutoff in set_snpnum:
        win_len = [num_cutoff] + list(windows[num_cutoff][x]['end'] - windows[num_cutoff][x]['start'] + 1 for x in windows[num_cutoff])
        f.write('\t'.join(list(str(x) for x in win_len)) + '\n')