# Calculate synonymous pairwise window (F<sub>ST</sub>) between populations
----
This program calculates pairwise population differentiation of allele frequenccies (F<sub>ST</sub>) across genome-wide 4-fold synonymous SNPs. 

## Environment
Follow instructions in [main documentation](../../README.md) to set up the conda environmnet '[WGS_analysis](../../conda_environment/WGS_analysis.yml)'.

## Code by step
1. Generate vcf files and mpileup files of 4-fold synonymous SNPs and analyzed sites from previously generated vcf and mpileup files: configure file paths and parameters in [run_extract_4_syn_sites.sh](code/run_extract_4_syn_sites.sh) and run it in bash. 
2. Calculate Fst from files generated using the same scripts for [genome-wide F<sub>ST</sub>](../Fst/): configure file paths and parameters in [run_calc_syn_Fst_reynolds.sh](code/run_calc_syn_Fst_reynolds.sh) and run it in bash.

In most cases, you do not have to modify the core python scripts which [filter for analyzed sites within window](code/filter_mpileup_by_cov_window.py) and [calculates F<sub>ST</sub> within each window](code/calc_Fst_reynolds_ungapped_efs.py), as well as the [paralleling shell script](code/run_calc_syn_Fst_reynolds.sh), which parallelly calculate within-window F<sub>ST</sub> and then merged results into genome-wide F<sub>ST</sub>.

Similarly, you do not have to modify the [core python script](code/extract_4_syn_sites.py) which extract 4-fold synonymous sites, and the [paralleling shell script](code/extract_4_syn_sites.sh), which parallely extract sites within each chromosome/contig.

## Input & output
Run the following commands to see input and output files for each step:
```
python extract_4_syn_sites.py -h
python calc_Fst_reynolds_ungapped_efs.py -h
```

You will need the output from [genomic annotation](../../genomic_annotation/) to extract 4-fold synonymous sites.