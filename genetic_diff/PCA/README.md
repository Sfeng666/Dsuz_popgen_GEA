#  Perform PCA on minor allele frequencies across populations

----
This program calculates biallelic minor allele frequency, and perform PCA to look at genetic differentiation among populations samples.

## Environment
Follow instructions in [main documentation](../../README.md) to set up the conda environmnet '[WGS_analysis](../../conda_environment/WGS_analysis.yml)'.

## Code by steps
1. Calculate a matrix of minor allele frequencies across SNPs and populations: configure file paths and parameters in this [shell script](code/run_calc_maf.sh) and run it in bash. 

    In most cases, you do not have to modify the [core python script](code/calc_maf_biallelic.py), which calculates MAF for SNPs within each contig, as well as the [paralleling shell script](code/calc_maf.sh), which parallelly calculates within-contig MAF matrices and then merged results into one matrix.

2. Perform PCA and plot samples on a 3D PCA space: configure file paths and parameters in this [shell script](code/run_PCA.sh) and run it in bash. 

    You could further modify the [core R script](code/plot_PCA.R) to change sample names and color options.

## Input & output
Check variables in this [shell script](code/run_calc_maf.sh) to see input and output files for calculating MAF.

Check variables in this [shell script](code/run_PCA.sh) and [R script](code/plot_PCA.R) to see input and output files for PCA.

## Notes
1. We calculated a matrix of allele frequencies across SNPs and population samples for PCA. Use any one of the two most frequent allele across all samples (sorted by the sum allele frequency from all samples, instead of allele count that adds weight to each voting sample with uneven depth). Choosing either Ref or alt doesn't make a difference as long as it's consistently chosen among samples.