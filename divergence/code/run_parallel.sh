#!/bin/bash

# run the job
bash parallel_by_block.sh \
-w /raid10/siyuan/marula_storage/jobs/suzukii_WGS/reannotate/calc_divergence/corrected_N_sites \
-a /raid10/siyuan/marula_storage/jobs/suzukii_WGS/reannotate/for_reannotate.fa \
-m /home/siyuan/reference/WT3-2.0/multi_align/suz_3.0_bia_1.0_mel_6.03/map_modified \
-s /raid10/siyuan/marula_storage/jobs/suzukii_WGS/reannotate/calc_divergence/corrected_N_sites/count_divergent_sites.py \
-r /home/siyuan/reference/WT3-2.0/refseq/annotation/GCF_013340165.1_LBDM_Dsuz_2.1.pri_assembly_report.txt \
-t /home/siyuan/reference/WT3-2.0/multi_align/suz_3.0_bia_1.0_mel_6.03/suz \
-n /home/siyuan/reference/WT3-2.0/multi_align/suz_3.0_bia_1.0_mel_6.03/bia \
-b /opt/miniconda3/envs/play_with_genome/bin/bedtools \
-d /opt/miniconda3/ \
-e WGS_analysis \
> run.log 2>&1