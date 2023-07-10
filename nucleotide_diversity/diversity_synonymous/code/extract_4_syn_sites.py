import os
from optparse import OptionParser
import gzip

### help and usage ### 
usage = "usage: %prog [options] args"
description = '''Function: extract 4-fold synonymous SNPs and filtered sites'''
version = '%prog 07.26.2021'
parser = OptionParser(usage=usage,version=version, description = description)
parser.add_option("--in_snp",
                    action="store",
                    dest = "in_snp",
                    help = "input path of filtered SNPs from a specified contig in vcf format (.vcf)",
                    metavar = "PATH")
parser.add_option("--in_mpileup",
                    action="store",
                    dest = "in_mpileup",
                    help = "input path of filtered sites from a specified contig in mpileup format (.gz)",
                    metavar = "PATH")
parser.add_option("--in_anno_seq",
                    action = "store",
                    dest = "in_anno_seq",
                    help = "input path of reannotation of the whole genome in fasta format (.fa)",
                    metavar = "PATH")
parser.add_option("--out_snp_syn",
                    action = "store",
                    dest = "out_snp_syn",
                    help = "output path of extracted 4-fold synonymous snps in vcf format (.vcf)",
                    metavar = "PATH")
parser.add_option("--out_mpileup_syn",
                    action = "store",
                    dest = "out_mpileup_syn",
                    help = "output path of extracted filtered sites in mpileup format",
                    metavar = "PATH")
parser.add_option("--temp",
                    action = "store",
                    dest = "temp",
                    help = "output directory of temporary intermediate files",
                    metavar = "PATH")
parser.add_option("--bedtools",
                    action = "store",
                    dest = "bedtools",
                    help = "path of bedtools",
                    metavar = 'PATH')                    
(options,args) = parser.parse_args()

### introduced variables ###
in_snp = options.in_snp
in_mpileup = options.in_mpileup
in_anno_seq = options.in_anno_seq
out_snp_syn = options.out_snp_syn
out_mpileup_syn = options.out_mpileup_syn
temp = options.temp
bedtools=options.bedtools

### 1. Get the annotation of SNPs and filtered sites
os.system('''awk 'BEGIN{FS="\\t"} $0 !~ /#/ {printf "%s\\t%s\\t%s\\t.\\t0\\t+\\n", $1, $2 - 1, $2}' ''' + in_snp + ' > ' + temp + '/' + in_snp.split('/')[-1] + '.bed')
os.system('{1} getfasta -s -tab -fi {2} -bed {3}/{0}.bed -fo {3}/{0}.txt'.format(in_snp.split('/')[-1], bedtools, in_anno_seq, temp))

os.system('''zcat '''+ in_mpileup + ''' | awk 'BEGIN{FS="\\t"} $0 !~ /#/ {printf "%s\\t%s\\t%s\\t.\\t0\\t+\\n", $1, $2 - 1, $2}' > ''' + temp + '/' + in_mpileup.split('/')[-1].replace('.gz', '') + '.bed')
os.system('{1} getfasta -s -tab -fi {2} -bed {3}/{0}.bed -fo {3}/{0}.txt'.format(in_mpileup.split('/')[-1].replace('.gz', ''), bedtools, in_anno_seq, temp))


### 2. Extract SNPs and filtered sites by degeneracy
with open(in_snp, 'r') as f, open('{1}/{0}.txt'.format(in_snp.split('/')[-1], temp), 'r') as f_anno, open(out_snp_syn, 'w') as f_syn:
    for line in f:
        if line.startswith('#'):
            f_syn.write(line)
            if line.startswith('#CHROM'):
                break

    for line1, line2 in zip(f, f_anno):
        anno = line2.strip().split('\t')[1]
        if anno == 'S':
            f_syn.write(line1)

with gzip.open(in_mpileup, 'rt') as f, open('{1}/{0}.txt'.format(in_mpileup.split('/')[-1].replace('.gz', ''), temp), 'r') as f_anno, open(out_mpileup_syn, 'w') as f_syn:
    for line1, line2 in zip(f, f_anno):
        anno = line2.strip().split('\t')[1]
        if anno == 'S':
            f_syn.write(line1)

os.system('rm {1}/{0}.bed {1}/{0}.txt'.format(in_snp.split('/')[-1], temp))
os.system('rm {1}/{0}.bed {1}/{0}.txt'.format(in_mpileup.split('/')[-1].replace('.gz', ''), temp))