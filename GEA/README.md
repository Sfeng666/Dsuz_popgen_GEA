# Perform genome-environment association (GEA) analysis

----
This program performs association analysis between environmental and genetic differentiation among populations, in order to detect local adaptation linked to a given environmental variable. To reduce the running time overall, SNPs will be subsampled along chromosomes and run on CHTC in parallel.

## Update
These test runs are to test the sensitivity of BayeScEnv to the prior paramater p for locus-specific model against local adaptation model, at 0.5 and 0.8. The parameter was set at 0.2 previously, generating an inplausibly large amount of outliers.

## Code by steps
1. Prepare input files:
    * Environmental data: environmental information in the form of ‘environmental differentiation’, which is computed as a contrast to a reference (usually the average environment, but not obligatory). For each environmental variable, the input is one line that has each column as the distance of each sample (cautious about the order!).
        * Open the R Markdown script [extract_env_data.Rmd](code/extract_env_data.Rmd) in Rstudio and run code blocks step by step.
    * Genotype data: allele count of SNPs that are subsampled from auto or X-linked chromosomes. 
        * Configure file paths in the shell script [run_prep_gt_data.sh](code/run_prep_gt_data.sh) and run it in bash. It will call the core python script [prep_gt_data_split.py](code/prep_gt_data_split.py) to prepare genotype data, which you typically don't need to modidy.
2. Run BayScEnv for subsamples of SNPs with each environmental variable parallelly on CHTC. This workflow was built on [DAGMan (Directed Acyclic Graph Manager)](https://htcondor.org/dagman/dagman.html), and is primarily designed to run through the [HTCondor](https://htcondor.org/htcondor/overview/) job scheduler (set up to run on UW-Madison's CHTC). 
    * First generate DAGMan job submission file for each environmental variable using the python script [prepare_dag_bayescenv_split.py](code/prepare_dag_bayescenv_split.py)
    * Then submit all job submission files using the bash script [submit_dag.sh](code/submit_dag.sh)
    * You typically don't have to modify the core shell script [run_bayescenv_chtc.sh](code/run_bayescenv_chtc.sh) and configure file [run_bayescenv_chtc.sub](code/run_bayescenv_chtc.sub), unless you want to chaneg job submission parameters or run the program on a different job platform.
3. Merge results from split BayeScEnv jobs and recalculate q value for each site across all contigs including autosomal and X-linked contigs. 
    * Run the python script [merge_split_bayescenv.py](code/merge_split_bayescenv.py)
4. Pare down closely linked candidate SNPs by maintaining the site with the lowest q within each 20 kb genomic window, and identify the closest gene in each direction within a 200-exon flanking region that overlapped with each candidate SNP.
    * First generate perl script for each environmental variable from a template perl script [GO_SNP_parallel_linkedGenes_modifiedsuz.pl](code/GO_SNP_parallel_linkedGenes_modifiedsuz.pl)  by running the 

## Notes
1. To make this GEA analysis computationally feasible with our large SNP set, while still analyzing all qualifying SNPs, we applied a split-run strategy: we subsampled SNPs across concatenated sequences of contigs within the autosomes and the X chromosome separately, and then ran subsamples with BayesScEnv in parallel.  Since the null model of population structure is estimated separately in each run, we subsampled non-adjacent SNPs at a fixed interval to limit locus-specific biases in that estimation, where the length of the interval between jointly analyzed SNPs was equal to the total number of subsamples. The subsampling procedure:
    * sample size = ~10000 SNPs
    * subsample with a fixed interval that will generate subsamples of the sample size evenly across the genome
    * the actual sample size may be slightly lower than the specialized one, but will be the same for most subsamples
    * separately subsample from concatenated X and autosomes.
2. Since GEA was run on split windows of SNPs, the original q-value (qval_g) was also computed locally. To get q-values for genome-wide screening, we recalculate q-value based on the ranking of sites across genome in Step 3.