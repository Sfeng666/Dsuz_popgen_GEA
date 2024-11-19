# -*- coding: utf-8 -*-
################################################################################
##  *
##  *  Function: annotate every site of a genome by exclusive categories of interest]
##  *  Author: Siyuan Feng
##  *  Mail: siyuanfeng.bioinfo@gmail.com
##  *  Version: 11.19.2024
##  *
################################################################################

# import timeit
import os
import re
import gzip
from optparse import OptionParser
# starttime = timeit.default_timer()

### help and usage ### 
usage = "usage: %prog [options] args"
description = '''Function: to annotate each site of a genome by nine pre-defined categories'''
version = '%prog 11.19.2024'
parser = OptionParser(usage=usage,version=version, description = description)
parser.add_option("--region",
                    action="store",
                    dest = "region",
                    help = "RefSeq Accession of the chromosome",
                    metavar = "ARG") 
parser.add_option("--in_anno",
                    action="store",
                    dest = "in_anno",
                    help = "input path of genome annotation in compressed gff3 format (.gff.gz)",
                    metavar = "PATH") 
parser.add_option("--in_genome_seq",
                    action = "store",
                    dest = "in_genome_seq",
                    help = "input path of genome sequence in fasta format (.fasta/.fna)",
                    metavar = "PATH")
parser.add_option("--in_condon_tab",
                    action = "store",
                    dest = "in_condon_tab",
                    help = "input path of a condon table in tab-delimited txt format (.txt)",
                    metavar = "PATH")
parser.add_option("--out_coord",
                    action="store",
                    dest = "out_coord",
                    help = "output path of coordinates of sites on both strands by categories in bed format (.bed)",
                    metavar = "PATH")
parser.add_option("--out_seq_for",
                    action = "store",
                    dest = "out_seq_for",
                    help = "output path of reannotated forward genome sequence of sites by categories in fasta format (.fasta)",
                    metavar = "PATH")
parser.add_option("--out_seq_rev",
                    action = "store",
                    dest = "out_seq_rev",
                    help = "output path of reannotated reverse genome sequence of sites by categories in fasta format (.fasta)",
                    metavar = "PATH")
parser.add_option("--out_report",
                    action = "store",
                    dest = "out_report",
                    help = "output path of reannotation report",
                    metavar = 'PATH')
parser.add_option("--bedtools",
                    action = "store",
                    dest = "bedtools",
                    help = "path of bedtools",
                    metavar = 'PATH')                                                                                                 
(options,args) = parser.parse_args()

### introduced variables ###
region = options.region
in_anno = options.in_anno
in_genome_seq = options.in_genome_seq
in_condon_tab = options.in_condon_tab
out_coord = options.out_coord
out_seq_for = options.out_seq_for
out_seq_rev = options.out_seq_rev
out_report = options.out_report
bedtools = options.bedtools

### in-script variables ###
annotation = {}
strands = ['+', '-']
coordinate = {std: {} for std in strands}
coordinate_CDS = {std: {} for std in strands}
prior = {'A': 9, 'W': 8, 'H': 7, 'S': 6, 'CDS': 6, '3': 5, '5': 4, 'R': 3, 'I': 2, 'G': 1}
prior_name = {'A': 'non-synonymous', 'W': '2-fold_synonymous', 'H': '3-fold_synonymous', 'S': '4-fold_synonymous', '3': "3'UTR", '5': "5'UTR", 'R': 'RNA-coding_regions', 'I': 'intron', 'G': 'intergenic'}

### functions ###

## determine the highest-priority feature among sites with overlap feature
def set_prior(list):
    return sorted(list, key=lambda x: prior[x], reverse=True)[0]

## determine the highest-priority feature among blocks with overlap feature when building the dictionary of coordinates
def set_prior_block(mark, dic, key):
    if key in dic:
        feature = set_prior([dic[key], mark])
    else:
        feature = mark
    return feature

## cut single strings of genome sequence into multiple lines at the length of 80
def cut_seq(seq, seq_len):
    lines = re.findall('.{' + str(seq_len) + '}', seq)
    lines.append(seq[(len(lines)*seq_len):])
    return lines

## function to replace a character at a specified location in a string
def sub(string, index, char):
    string = list(string)
    string[index] = char
    return ''.join(string)

## extract genomic sequence for concatenated CDS(s) of a gene and annotate the sequence as synonymous/non-synonymous
def extract_seq(path_bedtools, path_genome, CDS_dic, CDSs):
    # extract CDS sequences
    seq_index = {}
    for CDS in CDSs:
        coord = CDS_dic[CDS]
        chr = coord[0]
        start = str(int(coord[1]) - 1)
        end = coord[2]
        strand = coord[3]
        if seq_index == {}:
            seq_index[CDS] = [0, int(end) - int(start)]
        else:
            seq_index[CDS] = [seq_index[CDS - 1][1], seq_index[CDS - 1][1] + int(end) - int(start)]

        os.system('printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" {0} {1} {2} CDS 0 {3} >> {0}.bed'.format(chr, start, end, strand, path_bedtools, path_genome))
    
    os.system('{1} getfasta -s -tab -fi {2} -bed {0}.bed -fo {0}.txt'.format(chr, path_bedtools, path_genome))

    all_seq = ''
    with open('{0}.txt'.format(chr), 'r') as f:
        for line in f:
            all_seq += line.strip().split('\t')[1].upper()
    os.system('rm {0}.bed {0}.txt'.format(chr))
    # print(len(all_seq))

    # annotate CDS sequences
    phase_fstCDS = int(CDS_dic[0][4])
    all_anno_seq = 'A'*phase_fstCDS
    for codon in cut_seq(all_seq[phase_fstCDS:], 3):
        if len(codon) == 3:
            if 'N' in codon:
                all_anno_seq += 'AAA'
            else:           
                fold = {1: 'A', 2: 'W', 3: 'H', 4: 'S'}
                AA = codons[codon]                        
                for i in range(0, 3):
                    synonym = list(codons[sub(codon, i, base)] for base in bases).count(AA)
                    all_anno_seq += fold[synonym]
        else:
            all_anno_seq += 'A'*len(codon)
    return seq_index, all_anno_seq

### 1. Read in the annotation file from RefSeq gff
with gzip.open(in_anno, 'rt') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            chr = line[0]
            feature = line[2]
            start = line[3]
            end = line[4]
            strand = line[6]
            attr = line[8]
            phase = line[7]

            if chr == region:
                if feature == 'region':
                    annotation.setdefault(chr, dict(chr=[start, end]))
                    for std in strands:
                        coordinate[std].setdefault(chr, {})
                        coordinate_CDS[std].setdefault(chr, {})


                elif feature == 'gene' or feature == 'pseudogene':
                    gene = re.search('ID=(.*?);', attr).group(1)
                    try:
                        gene_biotype = re.search('gene_biotype=(.*?);', attr).group(1)
                    except:
                        gene_biotype = re.search('gene_biotype=(.*)', attr).group(1)
                    annotation[chr].setdefault(gene, dict(gene=[chr, start, end, strand, gene_biotype]))

                elif feature == 'mRNA' or feature == 'transcript' or feature in ["guide_RNA", "ncRNA", "miRNA", "lnc_RNA", "rRNA", "snRNA", "snoRNA", "SRP_RNA", "RNase_MRP_RNA", "RNase_P_RNA", "tRNA", "antisense_RNA", "primary_transcript"]:
                    RNA = re.search('ID=(.*?);', attr).group(1)
                    parent = re.search('Parent=(.*?);', attr).group(1)

                    if feature == 'mRNA':   # nest 'mRNA' under 'gene', and create empty lists of 'exon' and 'CDS' under 'mRNA'
                        annotation[chr][parent].setdefault(RNA, dict(exon = {}, CDS = {}))
                    else:   # nest other types of transcripts under 'gene', and create empty lists of 'exon' under the transcript
                        if feature != 'miRNA':  # the hierachy of miRNA can be nested under 'primary_transcript' in some cases, where there 'parent' is a 'primary_transcript' instead of a 'gene'. Since and the 'parent' transcript is not defined in annotation[chr], error would occur.
                            annotation[chr][parent].setdefault(RNA, dict(exon={}))

                elif feature == 'exon':
                    if gene_biotype == "pseudogene":
                        pass
                    else:
                        parent = re.search('Parent=(.*?);', attr).group(1)
                        if parent in annotation[chr][gene]: # added to avoid error caused by some exons that are nested under 'miRNA', which is further nested under 'primary_transcript', instead of under the level of a 'gene'
                            num = len(annotation[chr][gene][parent]['exon'])
                            annotation[chr][gene][parent]['exon'][num] = [chr, start, end, strand]
                            if not (gene_biotype == "protein_coding" and 'CDS' in annotation[chr][gene][parent]):
                                key = '#'.join([chr, start, end, strand])
                                coordinate[strand][chr][key] = set_prior_block('R', coordinate[strand][chr], key)
                    
                elif feature == 'CDS':
                    parent = re.search('Parent=(.*?);', attr).group(1)
                    num = len(annotation[chr][gene][parent]['CDS'])
                    annotation[chr][gene][parent]['CDS'][num] = [chr, start, end, strand, phase]

### 2. Add annotations of 5'UTR, 3'UTR, synonymous, non-synonymous sites, and intron
## read the codon table
bases = ['A', 'T', 'G', 'C']
codons = {}
with open(in_condon_tab, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            codons[line[0]] = line[2]

## annotate genomic elements
for chr in annotation:
    for gene in annotation[chr]:
        if gene != 'chr': 
            strand = annotation[chr][gene]['gene'][3]
            gene_biotype = annotation[chr][gene]['gene'][4]

            for RNA in annotation[chr][gene]:
                if RNA != 'gene' and RNA != 'exon' and RNA != 'CDS':

                    ## add annotations of intron
                    exons = sorted(annotation[chr][gene][RNA]['exon'], reverse=False)
                    exons_rev = sorted(annotation[chr][gene][RNA]['exon'], reverse=True)
                    if len(exons) > 1:
                        annotation[chr][gene][RNA].setdefault('intron', {})
                        if strand == "+":
                            exon_order = exons
                        else:
                            exon_order = exons_rev

                        for exon in exon_order:
                            start_intron = str(int(annotation[chr][gene][RNA]['exon'][exon][2]) + 1)
                            if exon_order.index(exon) != len(exon_order) - 1:
                                end_intron = str(int(annotation[chr][gene][RNA]['exon'][exon_order[exon_order.index(exon) + 1]][1]) - 1)
                            else:
                                break
                            num = len(annotation[chr][gene][RNA]['intron'])
                            annotation[chr][gene][RNA]['intron'][num] = [chr, start_intron, end_intron, strand]
                            key = '#'.join([chr, start_intron, end_intron, strand])
                            coordinate[strand][chr][key] = set_prior_block('I', coordinate[strand][chr], key)

                    if gene_biotype == 'protein_coding' and 'CDS' in annotation[chr][gene][RNA]:
                        ## add annotations of 5'UTR and 3'UTR
                        CDSs = sorted(annotation[chr][gene][RNA]['CDS'], reverse=False)
                        if strand == '+':
                            CDS_start = int(annotation[chr][gene][RNA]['CDS'][CDSs[0]][1])
                            CDS_end = int(annotation[chr][gene][RNA]['CDS'][CDSs[-1]][2])
                            exon_start = list(annotation[chr][gene][RNA]['exon'][exon][1] for exon in annotation[chr][gene][RNA]['exon'] 
                            if CDS_start >= int(annotation[chr][gene][RNA]['exon'][exon][1]) and CDS_start <= int(annotation[chr][gene][RNA]['exon'][exon][2]))[0]
                            exon_end = list(annotation[chr][gene][RNA]['exon'][exon][2] for exon in annotation[chr][gene][RNA]['exon'] 
                            if CDS_end >= int(annotation[chr][gene][RNA]['exon'][exon][1]) and CDS_end <= int(annotation[chr][gene][RNA]['exon'][exon][2]))[0]
                            annotation[chr][gene][RNA]['5UTR'] = [chr, exon_start, str(CDS_start - 1), strand]
                            annotation[chr][gene][RNA]['3UTR'] = [chr, str(CDS_end + 1), exon_end, strand]
                        else:
                            CDS_start = int(annotation[chr][gene][RNA]['CDS'][CDSs[0]][2])
                            CDS_end = int(annotation[chr][gene][RNA]['CDS'][CDSs[-1]][1])
                            exon_start = list(annotation[chr][gene][RNA]['exon'][exon][2] for exon in annotation[chr][gene][RNA]['exon']
                             if CDS_start >= int(annotation[chr][gene][RNA]['exon'][exon][1]) and CDS_start <= int(annotation[chr][gene][RNA]['exon'][exon][2]))[0]
                            exon_end = list(annotation[chr][gene][RNA]['exon'][exon][1] for exon in annotation[chr][gene][RNA]['exon']
                             if CDS_end >= int(annotation[chr][gene][RNA]['exon'][exon][1]) and CDS_end <= int(annotation[chr][gene][RNA]['exon'][exon][2]))[0]
                            annotation[chr][gene][RNA]['5UTR'] = [chr, str(CDS_start + 1), exon_start, strand]
                            annotation[chr][gene][RNA]['3UTR'] = [chr, exon_end, str(CDS_end - 1), strand]

                        if int(annotation[chr][gene][RNA]['5UTR'][1]) <= int(annotation[chr][gene][RNA]['5UTR'][2]):
                            key = '#'.join(annotation[chr][gene][RNA]['5UTR'])
                            coordinate[strand][chr][key] = set_prior_block('5', coordinate[strand][chr], key)

                        if int(annotation[chr][gene][RNA]['3UTR'][1]) <= int(annotation[chr][gene][RNA]['3UTR'][2]):
                            key = '#'.join(annotation[chr][gene][RNA]['3UTR'])
                            coordinate[strand][chr][key] = set_prior_block('3', coordinate[strand][chr], key)

                        ## add annotations of synonymous and non-synonymous sites
                        seq_index, all_anno_seq = extract_seq(bedtools, in_genome_seq, annotation[chr][gene][RNA]['CDS'], CDSs)
                        # print(seq_index, len(all_anno_seq))                        
                        for CDS in CDSs:
                            coord = annotation[chr][gene][RNA]['CDS'][CDS]
                            if strand == '+':
                                anno_seq = all_anno_seq[seq_index[CDS][0]: seq_index[CDS][1]]
                            else:
                                anno_seq = all_anno_seq[seq_index[CDS][0]: seq_index[CDS][1]][::-1]
                            # print(len(anno_seq), int(coord[2]) - int(coord[1]) + 1, annotation[chr][gene][RNA]['CDS'][CDS][-1])
                            key = '#'.join(coord[:-1])
                            coordinate[strand][chr][key] = set_prior_block('CDS', coordinate[strand][chr], key)

                            if coordinate[strand][chr][key] == 'CDS':
                                coordinate_CDS[strand][chr].setdefault(key, {})

                                for i in range(0, len(anno_seq)):
                                    if i > 0:
                                        if anno_seq[i] != anno_seq[i - 1]:
                                            try:
                                                consec_seq = re.search('^.*[^{0}](.*)'.format(anno_seq[i - 1]), anno_seq[:i]).group(1)
                                            except:
                                                consec_seq = anno_seq[:i]
                                            start = anno_seq[:i + 1].rindex(consec_seq)
                                            end = start + len(consec_seq) - 1
                                            key_site = '#'.join([chr, str(int(coord[1]) + start), str(int(coord[1]) + end), coord[3]])
                                            coordinate_CDS[strand][chr][key][key_site] = anno_seq[i - 1]
                                            if i == len(anno_seq) - 1:
                                                key_site = '#'.join([chr, str(int(coord[1]) + i), str(int(coord[1]) + i), coord[3]])
                                                coordinate_CDS[strand][chr][key][key_site] = anno_seq[i]
                                        else:
                                            if i == len(anno_seq) - 1:
                                                try:
                                                    consec_seq = re.search('^.*[^{0}](.*)'.format(anno_seq[i]), anno_seq[:i + 1]).group(1)
                                                except:
                                                    consec_seq = anno_seq[:i + 1]
                                                start = anno_seq[:i + 1].rindex(consec_seq)
                                                end = start + len(consec_seq) - 1
                                                key_site = '#'.join([chr, str(int(coord[1]) + start), str(int(coord[1]) + end), coord[3]])
                                                coordinate_CDS[strand][chr][key][key_site] = anno_seq[i]
                                    else:
                                        if len(anno_seq) == 1:
                                            key_site = '#'.join([chr, coord[1], coord[1], coord[3]])
                                            coordinate_CDS[strand][chr][key][key_site] = anno_seq

### 3. Mark the genome with defined categories and output the results
with open(out_coord, 'w') as f_coord, open(out_seq_for, 'w') as f_seq_for, open(out_seq_rev, 'w') as f_seq_rev, open(out_report, 'a') as f_repo:
    ## screen through each base within each chromosome to determine classifications on sites including ambiguous ones
    class_count_both = {x: 0 for x in prior if x != 'CDS'}
    if not os.path.exists(out_report):
        f_repo.write('\t'.join(['Categories', 'Base count', 'Percent']) + '\n')

    for strand in strands:
        class_count = {x: 0 for x in prior if x != 'CDS'}
        for chr in annotation:
            chr_start = annotation[chr]['chr'][0]
            chr_end = annotation[chr]['chr'][1]
            anno_seq = ''

            for base in range(int(chr_start), int(chr_end) + 1):
                marks = list(coordinate[strand][chr][x] for x in coordinate[strand][chr] if base >= int(x.split('#')[1]) and base <= int(x.split('#')[2]))
                if not marks:
                    high_prior = 'G'
                else:
                    high_prior = set_prior(marks)
                    if high_prior == 'CDS':
                        marks_CDS = []
                        region_CDS = list(x for x in coordinate_CDS[strand][chr] if base >= int(x.split('#')[1]) and base <= int(x.split('#')[2]))
                        for CDS in region_CDS:
                            marks_CDS += list(coordinate_CDS[strand][chr][CDS][x] for x in coordinate_CDS[strand][chr][CDS] if base >= int(x.split('#')[1]) and base <= int(x.split('#')[2]))
                        high_prior = set_prior(marks_CDS)
                anno_seq += high_prior
                class_count[high_prior] += 1
                class_count_both[high_prior] += 1

                if base > 1:
                    if anno_seq[-1] != anno_seq[-2]:
                        try:
                            consec_seq = re.search('^.*[^{0}](.*)[^{0}]'.format(anno_seq[-2]), anno_seq).group(1)
                        except:
                            consec_seq = re.search('^(.*)[^{0}]'.format(anno_seq[-2]), anno_seq).group(1)
                        start = anno_seq.rindex(consec_seq)
                        end = start + len(consec_seq)
                        f_coord.write('\t'.join([chr, str(start), str(end), prior_name[anno_seq[-2]], '0', strand]) + '\n')
                        if base == int(chr_end):
                            f_coord.write('\t'.join([chr, str(base - 1), str(base), prior_name[anno_seq[-1]], '0', strand]) + '\n')
                
                    else:
                        if base == int(chr_end):
                            try:
                                consec_seq = re.search('^.*[^{0}](.*)'.format(anno_seq[-1]), anno_seq).group(1)
                            except:
                                consec_seq = anno_seq
                            start = anno_seq.rindex(consec_seq)
                            end = start + len(consec_seq)
                            f_coord.write('\t'.join([chr, str(start), str(end), prior_name[anno_seq[-1]], '0', strand]) + '\n')
            if strand == '+':
                f_seq_for.write('>{0}\n'.format(chr))
                for line in cut_seq(anno_seq, 80):
                    f_seq_for.write(line + '\n')
            else:
                f_seq_rev.write('>{0}\n'.format(chr))
                for line in cut_seq(anno_seq, 80):
                    f_seq_rev.write(line + '\n')

        # f_repo.write('\t'.join(['Whole genome', str(sum(class_count.values())), '100%']) + '\n')
        for category in class_count:
            f_repo.write('\t'.join([prior_name[category] + '(' + strand + ')', str(class_count[category]), format(class_count[category]*100/sum(class_count.values()), '.2f') + '%']) + '\n')

    for category in class_count_both:
        f_repo.write('\t'.join([prior_name[category] + '(both)', str(class_count_both[category]), format(class_count_both[category]*100/sum(class_count_both.values()), '.2f') + '%']) + '\n')
