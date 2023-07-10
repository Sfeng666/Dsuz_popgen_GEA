# Calculate genome-wide nucleotide diversity
----
This program calculates nucleotide diversity across genome-wide SNPs.

## Environment
Follow instructions in [main documentation](../../README.md) to set up the conda environmnet '[WGS_analysis](../../conda_environment/WGS_analysis.yml)'.

## How to run this program
Configure file paths and parameters in this [shell script](code/run_paratest.sh) and run it in bash. 

In most cases, you do not have to modify the [core python script](code/filter_mpileup_by_cov.py), which counts the total number of analyzed sites, as well as the [paralleling shell script](code/calc_diversity_same_sites_ts_updated_poolsnp_paratest.sh), which parallelly calculates within-contig accumulative heterozygosity and then merged results into genome-wide nucleotide diversity.

## Input & output
Check variables in this [shell script](code/calc_diversity_same_sites_ts_updated_poolsnp_paratest.sh) to see input and output files.

## Notes
1. We adopted an unbiased estimator of nucleotide diversity (θ ̂_Π) based on heterozygosity (Π), which has been optimized for pool-seq data (Ferretti et al. 2013). 