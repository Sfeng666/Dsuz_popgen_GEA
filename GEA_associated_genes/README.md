# Identify environment-associated genes
----
This program identify genes associated with environment-associated adaptive loci from GEA.

## Code by step
1. This program pare down closely linked candidate SNPs (from [GEA](../GEA/)) by maintaining the site with the lowest q within each 20 kb genomic window, and identify the closest gene in each direction within a 200-exon flanking region that overlapped with each candidate SNP.
    * First generate perl script for each environmental variable from a template perl script [GO_SNP_parallel_genomic_permutation_modifiedsuz.pl](code/GO_SNP_parallel_genomic_permutation_modifiedsuz.pl) by running [modify_pl_linkedGene.py](code/modify_pl_linkedGene.py)
    * Then run [submit_jobs.sh](code/submit_jobs.sh) to identify associated genes
2. Generate tables of full information about significantly associated loci (including association statistics) and genes.
    * run [make_table_snp_gene_stats.py](code/make_table_snp_gene_stats.py)