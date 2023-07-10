from optparse import OptionParser

### help and usage ### 
usage = "usage: %prog [options] args"
description = '''Function: calculate the total number of pairwise differences between sequences from two populations '''
version = '%prog 08.08.2021'
parser = OptionParser(usage=usage,version=version, description = description)
parser.add_option("--in_snp",
                    action="store",
                    dest = "in_snp",
                    help = "input path of filtered SNPs from a specified contig in vcf format (.vcf)",
                    metavar = "PATH")
parser.add_option("--names",
                    action="store",
                    dest = "names",
                    help = "a comma-delimited string of sample names",
                    metavar = "STR")
(options,args) = parser.parse_args()

### introduced variables ###
in_snp = options.in_snp
names = options.names

### 1. Sum up pairwise differences and store the value in a dictionary-based matrix
samples = names.split(',')
diff_matrix = {sample: [0]*len(samples) for sample in samples}
with open(in_snp, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            ref = line[3].upper()
            alt = line[4].upper()
            alleles = [ref] + alt.split(',')
            afs = []

            for info in line[9:]:
                gt = list(int(x) for x in info.split(':')[0].split('/'))
                ct = [int(info.split(':')[1])] + list(int(x) for x in info.split(':')[2].split(','))

                if len(set(gt)) == 1:
                    af = {alleles[gt[0]]: 1}

                else:
                    allele_ct = {alleles[x]: ct[gt.index(x)] for x in gt}
                    allele_2 = sorted(allele_ct, key = lambda x: allele_ct[x], reverse = True)[:2]
                    allele_ct_2 = {x: allele_ct[x] for x in allele_2}
                    dp = sum(allele_ct_2.values())
                    af = {x: allele_ct_2[x]/dp for x in allele_ct_2}
                afs.append(af)

            for sample in samples:
                af_sp = afs[samples.index(sample)]
                diff_vect = []
                for af in afs:
                    diff = 0
                    for ale1 in af_sp:
                        for ale2 in af:
                            if ale1 != ale2:
                                diff += af_sp[ale1]*af[ale2]
                    diff_vect.append(diff)
                diff_matrix[sample] = list(x + y for x,y in zip(diff_matrix[sample], diff_vect))

### 2. Print the total pairwise differences in the form of a matrix table
header = ['Populations'] + samples
print('\t'.join(header))
for sample in samples:
    print('\t'.join([sample] + list(str(x) for x in diff_matrix[sample])))