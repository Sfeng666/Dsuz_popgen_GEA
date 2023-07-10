# Calculate genome-wide synonymous nucleotide diversity
----
This program calculates nucleotide diversity across genome-wide 4-fold synonymous SNPs, which is a way to minimize relative effects of sequencing errors on calculate diversity levels, as compared to all sites.

## Environment
Follow instructions in [main documentation](../../README.md) to set up the conda environmnet '[WGS_analysis](../../conda_environment/WGS_analysis.yml)'.

## How to run this program
Configure file paths and parameters in this [shell script](code/run_paratest.sh) and run it in bash. 

In most cases, you do not have to modify the core python scripts, which [counts the total number of analyzed sites](code/filter_mpileup_by_cov.py) and [extract 4-fold synonymous sites within each chromosome/contig](code/calc_dxy_contig.py), as well as the [paralleling shell script](code/calc_diversity_synonymous.sh), which parallelly calculates within-contig accumulative heterozygosity and then merged results into genome-wide nucleotide diversity.

## Input & output
Check variables in this [shell script](code/run_paratest.sh) to see input and output files.

## Detailed steps
1. Generate vcf files and mpileup files of 4-fold synonymous SNPs and filtered sites from previously generated vcf and mpileup files;
2. Calculate nucleotide diversity from files generated.