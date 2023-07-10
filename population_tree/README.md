# Build population tree and infer admixture from allele frequency

----
This program builds maximum-likelihood tree and infer admixture events from allele frequencies of SNPs across populations using [TreeMix](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002967).

## Environment
Follow instructions in [main documentation](../README.md) to set up the conda environmnet '[wgs_suzukii](../conda_environment/wgs_suzukii.yml)'.

## Code by step
1. Run TreeMix with different combinations of parameters:

	configure file paths and parameters in this [shell script](code/run_ML_tree.sh) and run it in bash.

	In most cases, you do not have to modify [the core python script](code/prep_input_forTreeMix_efs.py), which prepare the input table of SNP allele frequency, and the [paralleling shell script](code/TreeMix.sh), which parallelly run TreeMix with different combinations of parameters.

2. Plot population tree and admixture graph from TreeMix results, plot the curve of the fraction of explained variance at different numbers of migration events, and find the point where where the curve first reaches plateau:

	configure file paths and parameters in this [shell script](code/run_plot_ML_tree.sh) and run it in bash. 

	Similarly, you do not have to modify [the core R plotting script](code/plot_ML_tree.R), unless you wish to change the range of migration events to plot.

## Determine parameters
A common problem of TreeMix analysis is that the software itself cannot determine the most likely number of migration events (`-m` parameter). Also itself cannot determine the value of `-k` parameter, which accounts for linkage disequilibrium. In addition, we are not sure how thresholding minor allele frequency for the input SNP data would affect the inference. Therefore, we run TreeMix with a combination of following parameters:

* A. multiple runs with different number of migration events (0-20);
* B. Autosomal & X-linked loci;
* C. Heuristically choose a K such that the median consecutive length of a window is around 500 bp;
* D. Different minor allele threshold (1, 10);

The hierachical structure of combinations: D(=2) x B(=2) x C(=500) x A(=21)

## Detailed steps
1. Select a minor allele count threshold (D);
2. For each contig, 
	* calculate effective sample size at each site (instead of effective pool size, which is a component of calculating effective sample size, and is set as the actual pool size due to the lack of replicate population samples) based on [Ferretti et al, 2013](10.1111/mec.12522);
	* calculate window length at each level of the number of SNPs contained in a single window (1 - 500 SNPs);
	* adding the major and minor allele frequency of filtered SNPs across populations to an across-contig input SNP data file for TreeMix (X-linked or Autosomal-linked, B);
3. Across contigs (X-linked and Autosomal-linked combined), determine the value of K such that the median consecutive length of a window is closest to 500 bp (C);
4. Run TreeMix with the set of different number of migration events (A);
5. Plot ML trees and residual matrixes.
6. Once parameters B, C, D are confirmed, estimated the fraction of the variance in relatedness between populations that is accounted for by each model ($ f $) with a given number of migration events, and choose the models that have just reached the plateu of $ f $ to present.

## Notes
1. The required input SNP AF file of TreeMix doesn't contain contig information, which assumes that the order of the SNPs in the file is the order of the SNPs in the genome. We therefore only include SNPs which are multiples of the selected number of SNPs within each window.
2. Since SNPs are incontinuous, we set the boundries of windows in the following way: the starting position of the first window of each contig is 1, and the following boundries of each window are midpoint between the current ending SNP and the next starting SNP.
3. TreeMix takes longer than expected to complete when calculating standing error and P-value for each migration edge, and the running time of those jackknife estimates increases in a non-linear (probably exponential) way as the number of input loci increases. 