from optparse import OptionParser
from math import factorial
import sympy
import os

### help and usage ### 
usage = "usage: %prog [options] args"
description = '''Function: calculate window Fst (Reynolds' modification to Wright's original formulation) between two populations '''
version = '%prog 09.05.2021'
parser = OptionParser(usage=usage,version=version, description = description)
parser.add_option("--in_snp",
                    action="store",
                    dest = "in_snp",
                    help = "input path of filtered SNPs from a specified chr in vcf format (.vcf)",
                    metavar = "PATH")
parser.add_option("--in_mpileup",
                    action="store",
                    dest = "in_mpileup",
                    help = "input path of mpileup file of a specified chr in mpileup format (.mpileup)",
                    metavar = "PATH")
parser.add_option("--in_indel_pos",
                    action="store",
                    dest = "in_indel_pos",
                    help = "input path of the file containing indel positions (.txt)",
                    metavar = "PATH")
parser.add_option("--in_sample_size",
                    action="store",
                    dest = "in_sample_size",
                    help = "input path of allelic (auto/X) sample sizes of all samples analyzed (.txt)",
                    metavar = "PATH")
parser.add_option("--out_window_stats",
                    action="store",
                    dest = "out_window_stats",
                    help = "output path of window stats of all chrs for all samples analyzed (.txt)",
                    metavar = "PATH")
parser.add_option("--script",
                    action="store",
                    dest = "script",
                    help = "path of the python script counting the number of sites within each window (.py)",
                    metavar = "PATH")
parser.add_option("--window_heter",
                    action="store",
                    dest = "window_heter",
                    help = "the threshold of accumulative average heterozygosity among all samples analyzed, that is used to define windows",
                    metavar = "INT")
parser.add_option("--min_count",
                    action="store",
                    dest = "min_count",
                    help = "the threshold of minimum minor allele count (to exclude low-count alleles that might be sequencing errors)",
                    metavar = "FLOAT")                    
parser.add_option("--min_cov",
                    action="store",
                    dest = "min_cov",
                    help = "the threshold of minimum coverage",
                    metavar = "FLOAT")
parser.add_option("--max_cov",
                    action="store",
                    dest = "max_cov",
                    help = "input file containing the maximum coverage threshold of each contig",
                    metavar = "PATH")
parser.add_option("--contig",
                    action="store",
                    dest = "contig",
                    help = "a string indicating the contig being processed",
                    metavar = "STR")
parser.add_option("--chr",
                    action="store",
                    dest = "chr",
                    help = "a string indicating whether the contig being processed is auto- or X-linked",
                    metavar = "STR")
parser.add_option("--names",
                    action="store",
                    dest = "names",
                    help = "a comma-delimited string of sample names",
                    metavar = "STR")
parser.add_option("--temp",
                    action = "store",
                    dest = "temp",
                    help = "output directory of temporary intermediate files",
                    metavar = "PATH")
(options,args) = parser.parse_args()

### introduced variables ###
in_snp = options.in_snp
in_mpileup = options.in_mpileup
in_indel_pos = options.in_indel_pos
in_sample_size = options.in_sample_size
out_window_stats = options.out_window_stats
script = options.script
window_heter = int(options.window_heter)
min_cov = options.min_cov
max_cov = options.max_cov
min_count = options.min_count
contig = options.contig
chr = options.chr
names = options.names
temp = options.temp

### def functions ###

## function to calculate effective sample size as the expectation of unique chromosomal draws
def expectation_sympy(nr, nc):
    return float(sum(j*factorial(nc)*sympy.functions.combinatorial.numbers.stirling(nr, j)/(factorial(nc - j)*nc**nr) for j in range(1, nc + 1)))

## function to calculate the Reynolds Fst
def Fst_reynolds(p1_afs, p2_afs, size1, size2):
    # p1_afs and p2_afs are dicitonaries that include allele frequency spectrum in each population
    # size 1 and size 2 are allelic samples size of each population (the number of auto/X chromosomes, instead of the number of individuals)
    # simply implementing the Reynolds's formulation (Reynolds et al, 1983), while the variable names is indicating the location of terms in the formulations.

    # the first term that is shared by both numerator (al) and denominator (al + bl)
    first_term = sum(list((p1_afs[x] - p2_afs[x])**2 for x in p1_afs))/2

    # the expression at the numerator of the second term that is also shared by both numerator (al) and denominator (al + bl)
    second_term_num2 = size1*(1 - sum(list(x**2 for x in p1_afs.values()))) + size2*(1 - sum(list(x**2 for x in p2_afs.values())))

    # implement al
    al_frac_num = ((size1 + size2)/2)*second_term_num2
    al_frac_denom = size1*size2*((size1 + size2)/2- 1)
    frac = al_frac_num/al_frac_denom
    al = first_term - frac
    if al < 0:
        al = 0
 
    # implement al + bl
    albl_frac_num = (size1*size2 - (size1 + size2)/2)*second_term_num2
    albl_frac = albl_frac_num/al_frac_denom
    albl = first_term + albl_frac

    # seperately return al and al + bl for Fst calculation (weighted average of sites within windows)
    return al, albl

### 1. read the input of sample size
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
        
### 2. partition windows based on summed heterozygosity and calculate pairwise Fst within windows
windows = {}
win, sum_heter, num_snp = 0, 0, 0
windows[win] = {'start': 1, 
'Fst_matrix_num': {sample: [0]*len(samples) for sample in samples},
'Fst_matrix_denom': {sample: [0]*len(samples) for sample in samples}}

samples = names.split(',')
expectation_dict = {}
start = 0
end = 0
with open(in_snp, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            num_snp += 1
            pos = int(line[1])
            ref = line[3].upper()
            alt = line[4].upper()
            alleles = [ref] + alt.split(',')
            afs = []
            dps = []

            # specify the ending boundry of each window at the midpoint of the end SNP and the next starting SNP
            if end == 1 and start == 1:
                sum_pos = windows[win]['end'] + pos
                if sum_pos & 1 == 0:
                    windows[win]['end'] = int(sum_pos/2)
                else:
                    windows[win]['end'] = int((sum_pos - 1)/2)
                win += 1
                windows.setdefault(win, {'start': windows[win - 1]['end'] + 1,
                'Fst_matrix_num': {sample: [0]*len(samples) for sample in samples},
                'Fst_matrix_denom': {sample: [0]*len(samples) for sample in samples}})
                end = 0
                start = 0

            # specify the starting boundry of each window at the midpoint of the starting SNP and the previous ending SNP + 1
            if start == 1:
                windows.setdefault(win, {'start': windows[win - 1]['end'] + 1,
                'Fst_matrix_num': {sample: [0]*len(samples) for sample in samples},
                'Fst_matrix_denom': {sample: [0]*len(samples) for sample in samples}})
                start = 0                
            
            # extract allele frequency of all alleles for all samples
            for info in line[9:]:
                gt = list(int(x) for x in info.split(':')[0].split('/'))
                ct = [int(info.split(':')[1])] + list(int(x) for x in info.split(':')[2].split(','))
                dp = int(info.split(':')[3])

                if len(set(gt)) == 1:
                    af = {alleles[gt[0]]: 1}
                    
                else:
                    allele_ct = {alleles[x]: ct[gt.index(x)] for x in gt}
                    af = {x: allele_ct[x]/dp for x in allele_ct}

                for ale in alleles:
                    if not ale in af:
                        af[ale] = 0

                afs.append(af)
                dps.append(dp)
                sum_heter += 1 - sum(list(x**2 for x in af.values()))

            # calculate pairwise Fst for single sites
            for sample1 in samples:
                af1 = afs[samples.index(sample1)]
                dp1 = dps[samples.index(sample1)]
                Fst_vect_num = []
                Fst_vect_denom = []
                for sample2 in samples:
                    af2 = afs[samples.index(sample2)]
                    dp2 = dps[samples.index(sample2)]
                    key1 = ','.join([str(dp1), str(sample_size[chr][sample1])])
                    key2 = ','.join([str(dp2), str(sample_size[chr][sample2])])
                    if key1 in expectation_dict:
                        size1 = expectation_dict[key1]
                    else:
                        size1 = expectation_sympy(dp1, sample_size[chr][sample1])
                        expectation_dict[key1] = size1
                    if key2 in expectation_dict:
                        size2 = expectation_dict[key2]
                    else:
                        size2 = expectation_sympy(dp2, sample_size[chr][sample2])
                        expectation_dict[key2] = size2                    
                    Fst_frac = Fst_reynolds(af1, af2, size1, size2)
                    al = Fst_frac[0]
                    albl = Fst_frac[1]
                    Fst_vect_num.append(al)
                    Fst_vect_denom.append(albl)
                windows[win]['Fst_matrix_num'][sample1] = list(x + y for x,y in zip(windows[win]['Fst_matrix_num'][sample1], Fst_vect_num))
                windows[win]['Fst_matrix_denom'][sample1] = list(x + y for x,y in zip(windows[win]['Fst_matrix_denom'][sample1], Fst_vect_denom))

            # partition windows based on accumulative heterozygosity vs threshold
            if sum_heter/len(samples) >= window_heter:
                windows[win]['end'] = pos
                windows[win]['num_snp'] = num_snp
                end = 1
                start = 1

                sum_heter = 0
                num_snp = 0

    if sum_heter/len(samples) < window_heter:
        if sum_heter/len(samples) >= window_heter*0.95:
            windows[win]['end'] = pos
            windows[win]['num_snp'] = num_snp
        else:
            del windows[win]

### 3. output the window positions as the input of a script counting the number of filtered sites within each window
if len(windows) > 0:
    window_pos = contig + '.pos'
    window_ct = contig + '.ct'
    with open(window_pos, 'w') as f:
        for win in windows:
            f.write('\t'.join([str(windows[win]['start']), str(windows[win]['end'])]) + '\n')

    ### 4. run the script to count the number of filtered sites within each window
    ct_script = contig + '.sh'
    indelfree_mpileup = contig + '.mpileup'
    with open(ct_script, 'w') as f:
    #     f.write('''awk 'BEGIN{FS="\\t"; sites=0} \
    # NR==FNR {D[$1$2]++; next} \
    # !($1$2 in D) && NR>FNR {print}' \
    # <(cat ''' + in_indel_pos + ''') \
    # <(zcat ''' + in_mpileup + ''') \
    # | parallel \
    # --tmpdir ''' + temp + ''' \
    # -k \
    # --pipe \
    # -j 100% \
    # --no-notice \
    # --cat python2.7 ''' + script + ''' \
    # --mpileup {} \
    # --min-cov ''' + min_cov + ''' \
    # --max-cov ''' + max_cov + ''' \
    # --min-count ''' + min_count + ''' \
    # --min-freq 0 \
    # --miss-frac 0.001 \
    # --names ''' + names + ''' \
    # --windows ''' + window_pos + ''' \
    # --windows_ct ''' + window_ct)
        f.write('''awk 'BEGIN{FS="\\t"} \
NR==FNR {D[$1$2]++; next} \
!($1$2 in D) && NR>FNR {print}' \
<(cat ''' + in_indel_pos + ''') \
<(zcat ''' + in_mpileup + ''') \
> ''' + indelfree_mpileup + '''
python2.7 ''' + script + ''' \
--mpileup ''' + indelfree_mpileup + ''' \
--min-cov ''' + min_cov + ''' \
--max-cov ''' + max_cov + ''' \
--min-count ''' + min_count + ''' \
--min-freq 0 \
--miss-frac 0.001 \
--names ''' + names + ''' \
--windows ''' + window_pos + ''' \
--windows_ct ''' + window_ct)

    os.system('bash ' + ct_script)
    os.system('rm ' + ct_script)
    os.system('rm ' + indelfree_mpileup)
    os.system('rm ' + window_pos)

    ### 5. output the stats of windows
    with open(out_window_stats, 'a') as f, open(window_ct, 'r') as f_ct_win:
        i = 0
        for line in f_ct_win:
            if i == 0:
                num_sites = line.strip()
            else:
                line = line.strip().split('\t')
                windows[int(line[0])]['num_sites'] = int(line[3])
                # header - 'Chr' 'Contig' 'Window' 'Start' 'End' 'Length' 'Number of filtered sites' 'Number of SNPs'
                f.write('\t'.join([chr, contig, line[0], line[1], line[2], str(int(line[2]) - int(line[1]) + 1), line[3], str(windows[int(line[0])]['num_snp'])]) + '\n')
            i += 1
    os.system('rm ' + window_ct)

    ### 6. Print the weighted (by the number of filtered sites) pairwise Fst in the form of a matrix table
    header = [num_sites] + samples
    print('\t'.join(header))
    for sample in samples:
        Fst_vect = [0]*len(samples)
        for win in windows:
            Fst_vect_add = list(x*windows[win]['num_sites']/y for x,y in zip(windows[win]['Fst_matrix_num'][sample], windows[win]['Fst_matrix_denom'][sample]))
            Fst_vect = list(x + y for x,y in zip(Fst_vect_add, Fst_vect))
        print('\t'.join([sample] + list(str(x) for x in Fst_vect)))