#!/bin/bash

# run the job
bash parallel_by_chr.sh \
-w /home/siyuan/jobs/suzukii_WGS/reannotate/in_parallel \
-s /home/siyuan/jobs/suzukii_WGS/reannotate/in_parallel/genome_reannotate.py \
-b /opt/miniconda3/envs/play_with_genome/bin/bedtools \
-a /home/siyuan/reference/WT3-2.0/refseq/annotation/GCF_013340165.1_LBDM_Dsuz_2.1.pri_genomic.gff.gz \
-g /home/siyuan/reference/WT3-2.0/refseq/genome_seq/GCF_013340165.1_LBDM_Dsuz_2.1.pri_genomic.fna \
-c /home/siyuan/jobs/suzukii_WGS/reannotate/codon.txt \
-d /opt/miniconda3/ \
-e play_with_genome \
> run.log 2>&1 &