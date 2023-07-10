# Calculate synonymous between-population sequence distance (D<sub>XY</sub>)

----
This program calculates pairwise between-population sequence distance (D<sub>XY</sub>) across genome-wide 4-fold synonymous SNPs.

## Environment
Follow instructions in [main documentation](../../README.md) to set up the conda environmnet '[WGS_analysis](../../conda_environment/WGS_analysis.yml)'.

## Code by step
1. Generate vcf files and mpileup files of 4-fold synonymous SNPs and analyzed sites from previously generated vcf and mpileup files: configure file paths and parameters in [run_extract_4_syn_sites.sh](code/run_extract_4_syn_sites.sh) and run it in bash (this step use the same code as synonymous F<sub>ST</sub>). 
2. Calculate D<sub>XY</sub> from files generated using the same scripts for [genome-wide D<sub>XY</sub>](../Dxy): configure file paths and parameters in [run_calc_syn_dxy.sh](code/run_calc_syn_dxy.sh) and run it in bash.

In most cases, you do not have to modify the core python scripts which [counts the total number of analyzed sites](code/filter_mpileup_by_cov_window.py) and [calculates D<sub>XY</sub>](code/calc_dxy_contig.py) within each chromosome/contig, as well as the [paralleling shell script](code/run_calc_syn_dxy.sh), which parallelly calculate within-contig D<sub>XY</sub> and then merged results into genome-wide synonymous D<sub>XY</sub>.

Similarly, you do not have to modify the [core python script](code/extract_4_syn_sites.py) which extract 4-fold synonymous sites, and the [paralleling shell script](code/extract_4_syn_sites.sh), which parallely extract sites within each chromosome/contig.

## Input & output
Check variables in this [shell script](code/run_calc_dxy.sh) to see input and output files.

## Notes
1. D<sub>XY</sub> is an absolute measure of population differentiation, and is independent of levels of within-population diversity.
2. Though the formula is applicable to tri- or quadro-allelic sites, the calculation is limited to the two most frequent alleles for consistency with previous calculations of nucleotide diversity.
3. The two most frequent alleles doesn't have to be the same among samples, e.g., genotype A/T for sample1 and A/G for sample2, where dxy between sample1 and sample2 shoudl be calculated as AF(1A) x AF(2G) + AF(1T) x AF(2A) + AF(1T) x AF(2G). 