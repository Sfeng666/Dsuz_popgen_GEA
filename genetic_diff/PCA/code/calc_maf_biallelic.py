from optparse import OptionParser

### help and usage ### 
usage = "usage: %prog [options] args"
description = '''Function: calculate biallelic minor allele frequency'''
version = '%prog 08.15.2021'
parser = OptionParser(usage=usage,version=version, description = description)
parser.add_option("--in_snp",
                    action="store",
                    dest = "in_snp",
                    help = "input path of filtered SNPs from a specified contig in vcf format (.vcf)",
                    metavar = "PATH")
(options,args) = parser.parse_args()

### introduced variables ###
in_snp = options.in_snp

### 1. calculate biallelic minor allele frequency
with open(in_snp, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            ref = line[3].upper()
            alt = line[4].upper()
            alleles = {x: 0 for x in [ref] + alt.split(',')}
            alleles_ord = [ref] + alt.split(',')
            ct_sp = {}

            # record read count of all alleles of all samples into a dictionary, and add up the allele frequency of each allele across samples.
            for i in range(len(line[9:])):
                info = line[9:][i]
                gt = list(int(x) for x in info.split(':')[0].split('/'))
                ct = [int(info.split(':')[1])] + list(int(x) for x in info.split(':')[2].split(','))
                dp = int(info.split(':')[3])

                if len(set(gt)) == 1:
                    ct_sp[i] = {alleles_ord[gt[0]]: max(ct)}

                else:
                    ct_sp[i] = {alleles_ord[x]: ct[gt.index(x)] for x in gt}
                
                for ale in alleles:
                    if ale in ct_sp[i]:
                        alleles[ale] += ct_sp[i][ale]/dp

            # sort alleles by the total allele frequency, and determine the minor allele
            bialleles = sorted(alleles, key=lambda x: alleles[x], reverse=True)[:2]

            # extract MAF of each sample and print MAFs 
            mafs = []
            for i in range(len(line[9:])):
                dp_bi = sum(list(ct_sp[i][ale] for ale in bialleles if ale in ct_sp[i]))
                if bialleles[1] in ct_sp[i]:
                    if dp_bi != 0:
                        maf = ct_sp[i][bialleles[1]]/dp_bi
                    else:
                        maf = 0
                else:
                    maf = 0
                mafs.append(str(maf))
            print('\t'.join(mafs))