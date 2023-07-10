import os
from optparse import OptionParser

### help and usage ### 
usage = "usage: %prog [options] args"
description = '''Function: calculate divergence for sites of each annotation category '''
version = '%prog 07.27.2021'
parser = OptionParser(usage=usage,version=version, description = description)
parser.add_option("--blk",
                    action="store",
                    dest = "blk",
                    help = "block id of the macthed region (interger from 1-)",
                    metavar = "ARG")
parser.add_option("--in_assembly_rep",
                    action="store",
                    dest = "in_assembly_rep",
                    help = "input path of the NCBI assembly report of the target species",
                    metavar = "PATH")
parser.add_option("--in_anno_seq",
                    action = "store",
                    dest = "in_anno_seq",
                    help = "input path of reannotation of the whole genome of the target species in fasta format (.fa)",
                    metavar = "PATH")
parser.add_option("--in_map",
                    action="store",
                    dest = "in_map",
                    help = "input path of the map file that includes 0-based coordinates aligned genomic blocks",
                    metavar = "PATH")
parser.add_option("--in_mfa_tgt",
                    action="store",
                    dest = "in_mfa_tgt",
                    help = "input path of multiple fasta alignment file of matched blocks in the target assembly (.mfa)",
                    metavar = "PATH")
parser.add_option("--in_mfa_ans",
                    action="store",
                    dest = "in_mfa_ans",
                    help = "input path of multiple fasta alignment file of matched blocks in the ancestral (outgroup) assembly (.mfa)",
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
blk = options.blk
in_assembly_rep = options.in_assembly_rep
in_anno_seq = options.in_anno_seq
in_map = options.in_map
in_mfa_tgt = options.in_mfa_tgt
in_mfa_ans = options.in_mfa_ans
temp = options.temp
bedtools=options.bedtools

### in-script variables ###
ctg_name = {}
map_coord = {}
ans_sites = {}

### 1. Build a dictionary for contig names transferrings
with open(in_assembly_rep, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            ctg_name[line[0]] = line[6]

### 2. Build a dictionary for coordinates of matched blocks
with open(in_map, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        map_coord[line[0]] = [line[1], line[2], line[3], line[4]]

### 3. read the multiple fasta alignment file of this alignment, 
### and map bases of the outgroup species onto the genome of the target species.
pos = int(map_coord[blk][1])
match_seq = ''
with open(in_mfa_tgt, 'r') as f_tgt, open(in_mfa_ans, 'r') as f_ans, open('{0}/{1}.bed'.format(temp, blk), 'w') as f_bed:
    for line1, line2 in zip(f_tgt, f_ans):
        if not line1.startswith('>'):
            line1 = line1.strip()
            line2 = line2.strip() 
            for base1, base2 in zip(line1, line2):
                if base1 != '-':
                    pos += 1
                    if base2 != '-':
                        if base1 != 'N' and base2 != 'N':
                            f_bed.write('\t'.join([ctg_name[map_coord[blk][0]], str(pos - 1), str(pos), '.', '0', '+']) + '\n')
                            if base1 != base2:
                                match_seq += 's'
                            else:
                                match_seq += 'c'

### 4. Get the annotation of all matched sites
os.system('{0} getfasta -s -tab -fi {1} -bed {2}/{3}.bed -fo {2}/{3}.txt'.format(bedtools, in_anno_seq, temp, blk))

### 5. count the allele difference for each annotation category
anno_name = {'A': 'non-synonymous', 'W': '2-fold_synonymous', 'H': '3-fold_synonymous', 'S': '4-fold_synonymous', '3': "3'UTR", '5': "5'UTR", 'R': 'RNA-coding_regions', 'I': 'intron', 'G': 'intergenic'}
anno_divergence = {x: [0, 0] for x in anno_name}

with open('{0}/{1}.txt'.format(temp, blk), 'r') as f:
    for stat, line in zip(match_seq, f):
        anno = line.strip().split('\t')[1]
        anno_divergence[anno][1] += 1
        anno_divergence[anno][0] += stat.count('s')

ctgr_order = []
count_order = []
for ctgr in anno_divergence:
    ctgr_order.append(anno_name[ctgr])
    count_order += anno_divergence[ctgr]

print('\t'.join(ctgr_order))
print('\t'.join(list(str(x) for x in count_order)))




