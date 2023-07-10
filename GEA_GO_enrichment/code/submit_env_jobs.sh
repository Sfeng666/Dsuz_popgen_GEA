#!/bin/bash
cd /home/sfeng77/jobs/EAA/downstream/GSEA/random_SNP/scripts/output_gene_list/top_outliers/iterate/GSEA_top_500_genes/more_permutation_reps/results

# script to submit dags for each environmental variable
for dag in $(ls */split_*)
do
    condor_submit_dag $dag
done