# Identify autosomal and X-linked contigs
----
Identify autosomal and X-linked contigs based on correlation with average depth of contigs of known assignment.

## Environment
Follow instructions in [main documentation](../README.md) to set up the conda environmnet '[WGS_analysis](conda_environment/WGS_analysis.yml)'.

## Code by step
1. Calculate the average depth of each contig in each sample: [code](code/calc_mean_depth.sh)
2. Format the table of calculated average depth as well as previous assignment of another approach,  for step 3: [code](code/format_delim.sh)
3. Assign unplaced contigs to autosome/X-chromosome based on correlation with average depth of autosome/X-chromosome: [code](code/calc_cor_ambigous.Rmd)
4. (optional) Calculate the total length of unplaced contigs after chromosome assignment using different approach: [code](code/calc_length_amb_contigs.sh)