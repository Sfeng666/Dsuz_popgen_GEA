# Calculate and plot window nucleotide diversity across each chromosome

----
This program calculates amd plots nucleotide diversity in windows across each chromosome, in order to examine patterns of polymorphism across chromosome arms, and to indirectly reflect regions of different recombination rates, as well as to investigate whether our demographic inference might be skewed by large regions of low recombination.

## Environment
Follow instructions in [main documentation](../../README.md) to set up the conda environmnet '[WGS_analysis](../../conda_environment/WGS_analysis.yml)'.

## Code by step
1. Calculate window nucleotide diversity:

	configure file paths and parameters in this [shell script](code/run_calc_diversity_window.sh) and run it in bash.

	In most cases, you do not have to modify the core python scripts, which [count the number of analyzed sites within each window](code/win_div_by_filtered_sites.py) and [calculates window nucleotide diversity for windows within each contig](code/calc_diversity_window.py), as well as the [paralleling shell script](code/calc_diversity_window.sh), which parallelly run calculation for all contigs.

2. Plot window diversity across chromosome arms:

	configure file paths and parameters in this [shell script](code/run_plot_win_diversity.sh) and run it in bash. 

	Similarly, you do not have to modify [the core R plotting script](code/plot_win_diversity.R), unless you wish to change which populations to show in the plot.

## Detailed steps
1. Divide each contig into windows, where each window is a continuous genomic region that includes 125,000 analyzed sites;
2. Within each window, calculate nucleotide diversity;
3. Plot the distribution of window diversity of each chromosome arm for selected samples.

## Notes
1. The window is defined as a continuous genomic region that include 125000 filtered sites.
2. We cannot confirm the order of contigs along each chromosome arms, but it seems ok to simply concatenate the windows from contigs that could be mapped to each chromosome arm, other wise many of those contigs will have only several windows, which is less informative. Contigs of a chromosome arm are ordered by their total length.
3. The position of each window nucleotide diversity is not shown in the plot, because a continuous axis of genomic coordinates cannot be applied to concatenated contigs without knowing the order. Instead, the windows are only shown in the following order: (contig_1: window_1, window_2 ... window_n) (contig_2: window_1, window_2 ... window_n) ... (contig_n: window_1, window_2 ... window_n). 