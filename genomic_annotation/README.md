# Reannotate genomic features

Reannotate the whole genome by categories of genomic elements or site degeneracy, and generate annotations in both `.fasta` and `.bed` formats.
----

## Environment
Follow instructions in [main documentation](../README.md) to set up the conda environmnet '[play_with_genome](../conda_environment/play_with_genome.yml)'.

## How to run this program
Configure file paths in the [shell script](code/run_parallel.sh) and run it in bash. 

In most cases, you do not have to modify the [core python script](code/genome_reannotate.py), which runs genomic annotation for each chromosome/contig, and the [core shell script](code/parallel_by_chr.sh), which parallelly performs genomic annotation for all chromosome/contig and then merge them into one output.

## Input & output
Run the following command to see input and output files:
```
python genome_reannotate.py -h
``` 

## This program annotates each site of a genome by nine pre-defined categories
1. 2-fold synonymous sites (of CDS, protein-coding gene) - reannotated as *W* ;
2. 3-fold synonymous sites (of CDS, protein-coding gene) - reannotated as *H* ;
3. 4-fold synonymous sites (of CDS, protein-coding gene) - reannotated as *S*;
4. Non-synonymous sites (of CDS, protein-coding gene) - reannotated as *A*;
5. 5' UTR (of protein-coding gene) - reannotated as *5*;
6. 3' UTR (of protein-coding gene) - reannotated as *3*;
7. Intron (of protein-coding gene/non-proteincoding gene) - reannotated as *I*;
8. RNA-coding regions (non-protein-coding gene - exons) - reannotated as *R*;
9. Intergenic regions  - reannotated as *G*.

## How we deal with some tricky issues in this program
1. This script takes into account overlapping genes, which can cause ambiguouty in classifications of sites;
2. This script takes into account genes with multiple transcripts. It's worth noting that even protein-coding genes could have non-proteincoding transcripts. In addition to protein-coding genes, non-coding RNAs, such as lncRNA, miscRNA, etc, could also have multiple transcripts;
3. The way we deal with ovalapping genes and multiple-transcripts genes is to follow a priority of annotated categories in the following order:    
_Nonsynonymous > 2-fold synonymous > 3-fold synonymous > 4 fold synonymous >  3’ UTR > 5’ UTR > RNA-coding regions > Introns > Intergenic_;
4. This script takes into account partial annotations of CDS (5' or 3' partial). That means they lack annotations of 5' or 3' UTR, correspondingly. Also, this may lead to incomplete first codons;
5. This script considers misc_RNA (defined as "Transcribed untranslated RNA that no-one really knows what it does" by Ensemble) as RNA-coding regions (non-protein-coding), despite of it being annotated as "protein_coding" under the "gene_biotype" attribute in RefSeq annotations. It's worth noting that some of the genes of misc_RNA have multiple transcripts;
6. The START and STOP codons are included in the CDS. That is, if the locations of the start and stop codons are known, the first three base pairs of the CDS should correspond to the start codon and the last three correspond the stop codon;
7. Forward and reverse strands are annotated separately in `.fasta` format, which is consistent with the RefSeq annotation.
8. For better performance, we use [GNU parallel][parallel] to parallelly run jobs by each chromosome/contig.
[parallel]: https://www.gnu.org/software/parallel/

## How we deal with inperfect annotations
1. degeneracy of N bases? - assigned as non-synonymous
2. degeneracy of incomplete codons? - assign partial end as non-synonymous, and concatenate CDSs to make middle-end codons complete
3. degeneracy of CDS on negative strand? - seperate annotations on +/- strand