# Calculate genome-wide between-population sequence distance (D<sub>XY</sub>)

----
This program calculates pairwise between-population sequence distance (D<sub>XY</sub>) across genome-wide SNPs.

## Environment
Follow instructions in [main documentation](../../README.md) to set up the conda environmnet '[WGS_analysis](../../conda_environment/WGS_analysis.yml)'.

## How to run this program
Configure file paths and parameters in this [shell script](code/run_calc_dxy.sh) and run it in bash. 

In most cases, you do not have to modify the core python scripts, which [counts the total number of analyzed sites](code/filter_mpileup_by_cov.py) and [calculates Dsub>XY</sub> within each chromosome/contig](code/calc_dxy_contig.py), as well as the [paralleling shell script](code/calc_dxy.sh), which parallelly calculates within-contig between-population heterozygosity and then merged results into genome-wide D<sub>XY</sub>.

## Input & output
Check variables in this [shell script](code/run_calc_dxy.sh) to see input and output files.

## Notes
1. D<sub>XY</sub> is an absolute measure of population differentiation, and is independent of levels of within-population diversity.
2. Though the formula is applicable to tri- or quadro-allelic sites, the calculation is limited to the two most frequent alleles for consistency with previous calculations of nucleotide diversity.
3. The two most frequent alleles doesn't have to be the same among samples, e.g., genotype A/T for sample1 and A/G for sample2, where dxy between sample1 and sample2 shoudl be calculated as AF(1A) x AF(2G) + AF(1T) x AF(2A) + AF(1T) x AF(2G). 