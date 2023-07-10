import sys
from collections import defaultdict as d
import re
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun
#Modified by: Siyuan Feng

#########################################################   HELP   #########################################################################
usage="""python %prog \
      --mpileup data.mpileup \
      --min-cov 10 \
      --max-cov data.cov \
      --min-count 10 \
      --min-freq 0.01 \
      --mis-frac 0.1 \
      --base-quality-threshold 15 \
      --names Kib32,Tam10 \
      --coding 1.8 \
      > output.vcf"""
description = '''Function: originally used to call SNP, modified to count the number of filtered sites within each window '''
parser = OptionParser(usage=usage, description=description)
helptext="""
H E L P :
_________
"""
group=OptionGroup(parser,helptext)
#########################################################   parameters   #########################################################################

parser.add_option("--mpileup", dest="m", help="A mpileup file")
parser.add_option("--min-cov", dest="minc", help="The minimum coverage threshold: e.g. 10",default=10)
parser.add_option("--max-cov", dest="max", help="An input file with precomputed coverage thresholds")
parser.add_option("--min-count", dest="mint", help="The minimum number of counts of the alternative allele across all samples pooled",default=3)
parser.add_option("--min-freq", dest="minf", help="The minimum Frequency of the alternative allele across all samples pooled",default=0.01)
parser.add_option("--miss-frac", dest="mis", help="The minimum Frequency of the alternative allele across all samples pooled",default=0.1)
parser.add_option("--base-quality-threshold", dest="b", help="The Base-quality threshold for Qualities encoded in Sanger format (Illumina 1.8 format)",default=15)
parser.add_option("--names", dest="n", help="a comma separted list of thenames of all samples in the mpileup file")
parser.add_option("--coding", dest="c", help="the Illumina FASTQ quality coding",default=1.8)
parser.add_option("--win_len", dest="win_len", help="input paarameter of the numebr of filtered sites that should be contained in each window")
parser.add_option("--win_pos", dest="win_pos", help="output path of the file containing genomic positions of each window, as well as the number of analyzed sites within each")

parser.add_option_group(group)
(options, args) = parser.parse_args()


################################### functions ######################################

def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    import gzip
    if x=="-":
        y=sys.stdin
    elif x.endswith(".gz"):
        y=gzip.open(x,"r")
    else:
        y=open(x,"r")
    return y

def keywithmaxvalue(d):
    ''' This function resturns the key for the maximum value in a dictionary'''
    newhash=d(list)
    for k,v in d.items():
        newhash[v].append(k)
    return newhash[max(newhash.keys())]


def splitter(l, n):
    ''' This generator function returns equally sized cunks of an list'''
    #credit: Meric Lieberman, 2012
    i = 0
    chunk = l[:n]
    while chunk:
        yield chunk
        i += n
        chunk = l[i:i+n]

def extract_indel(l,sign):
    ''' This function returns an Indel from a sequence string in a pileup'''
    position = l.index(sign)
    numb =""
    i = 0
    while True:
        if l[position+1+i].isdigit():
            numb+= l[position+1+i]
            i+=1
        else:
            break

    seqlength = int(numb)
    sequence = l[position:position+i+1+seqlength]
    indel=sequence.replace(numb,"")

    return sequence,indel

################################## parameters ########################################
data=options.m
minimumcov=int(options.minc)
minimumcount=int(options.mint)
minimumfreq=float(options.minf)
missfrac=float(options.mis)
baseqthreshold=int(options.b)
phred=float(options.c)
win_len=int(options.win_len)
win_pos=options.win_pos

############################ calculate PHRED cutoff  #############################

# calculate correct PHRED score cutoff: ASCII-pc

if phred >=1.0 and phred <1.8:
    pc=64
else:
    pc=33

############################ get MAX coverage threshold  #############################
maximumcov=d(list)
for l in open(options.max,"r"):
    if l.startswith("#") or l.startswith("calculating"):
        continue
    k,v=l.split("\t")
    maximumcov[k]=map(int,v.split(","))
#print maximumcov

############################ get window positions  #############################
windows = {}
win, num_sites = 0, 0
windows[win] = {'start': 1, 
'sites': 0}

start = 0
end = 0
############################ parse MPILEUP ###########################################

# parse mpileup and store alternative alleles:
for line in load_data(data):    
    if len(line.split("\t"))<2:
        continue

    k = line[:-1].split('\t')
    CHR,POS,REF = k[:3]

    # specify the ending boundry of each window at the midpoint of the end SNP and the next starting SNP
    # specify the starting boundry of each window at the midpoint of the starting SNP and the previous ending SNP + 1
    if end == 1 and start == 1:
        sum_pos = windows[win]['end'] + int(POS)
        if sum_pos & 1 == 0:
            windows[win]['end'] = int(sum_pos/2)
        else:
            windows[win]['end'] = int((sum_pos - 1)/2)
        win += 1
        windows.setdefault(win, {'start': windows[win - 1]['end'] + 1,
        'sites': 0})
        end = 0
        start = 0

    ## only keep chromosomal arms with maximum coverage threshold
    if CHR not in maximumcov:
        #print CHR
        continue

    div = list(splitter(k,3))
    libraries=div[1:]
    # loop through libraries
    totalalleles=d(int)
    alleles=d(lambda:d(int))

    for j in range(len(libraries)):
        alleles[j]
        nuc = libraries[j][1]
        qualities = libraries[j][2]

        # test if seq-string is empty
        if nuc=="*":
            continue

        # find and remove read indices and mapping quality string
        nuc = re.sub(r'\^.',r'',nuc)
        nuc = nuc.replace('$','')

        # find and remove InDels
        while "+" in nuc or "-" in nuc:
            if "+" in nuc:
                insertion,ins=extract_indel(nuc,"+")
                nuc=nuc.replace(insertion,"")
            else:
                deletion,dele=extract_indel(nuc,"-")
                nuc=nuc.replace(deletion,"")


        # test for base quality threshold (if below: ignore nucleotide)
        #print len(nuc),len(qualities)
        nuc = "".join([nuc[x] for x in range(len(nuc)) if ord(qualities[x])-pc>=baseqthreshold])
        nuc = "".join([nuc[x] for x in range(len(nuc)) if nuc[x]!="*"])

        # read all alleles
        for i in range(len(nuc)):

            # ignore single nucleotide deletions
            if nuc[i]=="*":
                continue
            # count nucleotides similar to reference base
            if nuc[i] =="," or nuc[i] == ".":
                totalalleles[REF]+=1
                alleles[j][REF]+=1
                continue
            # count alternative nucleotides
            totalalleles[nuc[i].upper()]+=1
            alleles[j][nuc[i].upper()]+=1

    ## test if SNPs pass minimum count / minimum frequency threshold:
    for allele,counts in totalalleles.items():
        if counts<minimumcount or counts/float(sum(totalalleles.values()))<minimumfreq:
            del totalalleles[allele]

    ## create output for VCF
    ADP=sum(totalalleles.values())/len(libraries)
    ALT=[]
    ## set alternative allele order:
    for i in ["A","T","C","G"]:
        if i==REF:
            continue
        if i not in totalalleles:
            continue
        ALT.append(i)


    ## set ADP,NC,GT,AD and DP
    ADP=sum(totalalleles.values())/len(libraries)
    miss=0

    for j in range(len(libraries)):
        ## make empty entry if no allele counts for sample
        if j not in alleles:
            miss+=1
            continue

        alleleh = alleles[j]
        # remove alleles not counted in all samples
        for k,v in alleleh.items():
            if k != REF and k not in ALT:
                del alleleh[k]
        GT,AD,RD,FREQ,NC=[],[],0,[],0
        DP=sum(alleleh.values())

        ## test if mincoverage is still reached or if the maxcoverage is still exceeded when removing alleles that do not fullfill criteria; make empty entry if sample not fullfilling min/max coverage threshold
        if DP<minimumcov or DP>maximumcov[CHR][j]:
            miss+=1
            continue

    ## test if missing fraction of samples smaller than threshold:
    if miss/float(len(libraries))>missfrac:
        #print CHR,POS,"missing fraction",miss/float(len(libraries))
        continue

    ## if passed all filterings, add a analyzed site
    num_sites += 1

    # partition windows based on accumulative heterozygosity vs threshold
    if num_sites >= win_len:
        windows[win]['end'] = int(POS)
        windows[win]['sites'] = num_sites
        end = 1
        start = 1

        num_sites = 0

if num_sites < win_len:
    if num_sites >= win_len*0.95:
        windows[win]['end'] = int(POS)
        windows[win]['sites'] = num_sites
    else:
        del windows[win]


## output the position and the number of analyzed sites of each window
with open(win_pos, 'w') as f:
    for win in windows:
        f.write('\t'.join(list(str(x) for x in [win, windows[win]['start'], windows[win]['end'], windows[win]['sites']])) + '\n')