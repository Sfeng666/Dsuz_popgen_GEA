# Calculate divergence for sites of each annotation category  

----
This program calculates divergence between two species for sites of each annotated category.

## Environment
Follow instructions in [main documentation](../README.md) to set up the conda environmnet '[WGS_analysis](../conda_environment/WGS_analysis.yml)'.

## How to run this program
Configure file paths in the [shell script](code/run_parallel.sh) and run it in bash. 

In most cases, you do not have to modify the [core python script](code/count_divergent_sites.py), which count the number of differing reference alleles within each alignment block, and the [core shell script](code/parallel_by_block.sh), which parallelly count the number of different reference alleles within each alignment block, and calculate the proportion of differing reference alleles across all blocks for each site category.

## Input & output
Run the following command to see input and output files:
```
python count_divergent_sites.py -h
``` 

## Detailed steps
1. Get the coordinates of an alignment block in the target species (i.e. *Drosophila. suzukii*) from the 'map' file;
2. Read the multiple fasta alignment (.mfa) file of this alignment, and map bases of the outgroup species (i.e. *Drosophila. biarmipes*) onto the genome of the target species. This will generate a bed file of coordinates of matched sites on the target species, and a sequence of allele substitution/consistency (labeled as s/c);
3. Use the bed file to generate a tab-delimited table of site annotations;
4. Simultaneusly read the site annotation table and the sequence of allele substitution/consistency, while count the allele difference for each annotation category;
5. Calculate divergence as the proportion of differentiated sites for each annotation category;
6. (Optional) Plot divergence for each annotation category.

Step 1-4 is done by [count_divergent_sites.py](code/count_divergent_sites.py) for each alignment block, and step 5 is implemented by [parallel_by_block.sh](code/parallel_by_block.sh). Step 6 is implemented by [divergence_barplot.R](code/divergence_barplot.R).