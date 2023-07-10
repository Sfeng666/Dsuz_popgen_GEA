# GO enrichment of the top environment-associated genes

----
This program performs GSEA on top 500 genes via genomic permutation of outlier SNP positions (100,000,000 replicates), which accounts for the variability of gene length and the clustering of functionally related genes, as described in previous work (Pool et al. 2017)

## Code by steps
1. Run GSEA in 100000 parallel replicates on CHTC
    * First generate perl script for each environmental variable from a template perl script [GO_SNP_parallel_genomic_permutation_modifiedsuz.pl](code/GO_SNP_parallel_genomic_permutation_modifiedsuz.pl) by running [modify_pl.py](code/modify_pl.py)
    * Then generate DAGMan job submission file for each environmental variable by running [prepare_dag_GSEA.py](code/prepare_dag_GSEA.py)
    * Submit all job submission files by running the bash script [submit_dsubmit_env_jobsag.sh](code/submit_env_jobs.sh)
    * You typically don't have to modify the core shell script [run_GSEA_chtc.sh](code/run_GSEA_chtc.sh) and condor job submission file [run_GSEA_chtc.sub](code/run_GSEA_chtc.sub), unless you want to chaneg job submission parameters or run the program on a different job platform.
2. Merge permutation p values across GSEA permutation replicates
    * run [merge_rep_GSEA.py](code/merge_rep_GSEA.py)