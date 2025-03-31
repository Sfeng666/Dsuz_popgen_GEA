# Dsuz_popgen_GEA: population genomics and GEA analyses on *Drosophila suzukii*
----
Code for major analyses performed in population genomic analyses and genome-environment association (GEA) analyses on *Drosophila suzukii*.

Related publication:

[Feng S, DeGrey SP, Guédot C, Schoville SD, Pool JE. 2024. Genomic Diversity Illuminates the Environmental Adaptation of Drosophila suzukii. Genome Biology and Evolution 16:evae195.](https://doi.org/10.1093/gbe/evae195)

## Analyses list
1. [Identify autosomal and X-linked contigs](assign_contig/README.md)
2. [Annotate genomic features](genomic_annotation/README.md)
3. [Calculate sequence divergence between two species](divergence/README.md)
4. Estimate [genome-wide](nucleotide_diversity/diversity/README.md) and [synonymous](nucleotide_diversity/diversity_synonymous/README.md) nucleotide diversity, and [illustrate window diversity across chromosome arms](nucleotide_diversity/window_diversity/README.md)
5. Estimate [genome-wide](genetic_diff/Fst/README.md) and [synonymous](genetic_diff/Fst_synonymous/README.md) F<sub>ST</sub>
6. Estimate [genome-wide](genetic_diff/Dxy/README.md) and [synonymous](genetic_diff/Dxy_synonymous/README.md) D<sub>XY</sub>
7. [Perform PCA on allele frequencies](genetic_diff/PCA/README.md)
8. [Build population tree and infer admixture from allele frequency](population_tree/README.md)
9. [Perform genome-environment association (GEA) analysis](GEA/README.md)
10. [Identify environment-associated genes](GEA_associated_genes/README.md)
11. [Perform GO enrichment of the top environment-associated genes](GEA_GO_enrichment/README.md)

## Environment setup

To set up the environment for above analyses, you could use conda:
```
conda env create -n environemnt_name --file environemnt_name.yml
```
Environment needed for each analysis is indicated in documentation of each analysis. Environmental package lists (`environemnt_name.yml`) for cross-platform install could be found under [this directory](conda_environment).

If you have not installed conda, run the following command:
```
# download miniconda
curl -sL \
  "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" > \
  "Miniconda3.sh"
```
```  
# install miniconda
bash Miniconda3.sh
```
