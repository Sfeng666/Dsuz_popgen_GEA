# Calculate genome-wide pairwise window (F<sub>ST</sub>) between populations
----
This program calculates pairwise population differentiation of allele frequenccies (F<sub>ST</sub>) across genome-wide SNPs.

## Environment
Follow instructions in [main documentation](../../README.md) to set up the conda environmnet '[WGS_analysis](../../conda_environment/WGS_analysis.yml)'.

## How to run this program
Configure file paths and parameters in this [shell script](code/run_calc_Fst_reynolds.sh) and run it in bash. 

In most cases, you do not have to modify the core python scripts, which [filter for analyzed sites within window](code/filter_mpileup_by_cov_window.py) and [calculates F<sub>ST</sub> within each window](code/calc_Fst_reynolds_ungapped_efs.py), as well as the [paralleling shell script](code/calc_Fst_reynolds.sh), which parallelly calculate within-window F<sub>ST</sub> and then merged results into genome-wide F<sub>ST</sub>.

## Input & output
Run the following command to see input and output files:
```
python calc_Fst_reynolds_ungapped_efs.py -h
``` 

## Notes
1. The calculation is separately performed for autosomes and X-chromosome, to account for different allelic sample sizes, demogpraphic histories, and outcomes of natural selection.
2. The formulation being used here is a a modification (also called G<sub>ST</sub>) to Wright's original formulation (Reynolds 1983), also known as Reynold's estimator (theta) of the coancestry coefficient, which takes into account unequal samples sizes among populations, as well as sites with more than two alleles.
3. The formulation is applicable to tri- or quadro-allelic sites, but it also raises the issue that two populations might shares no allele but will have Fst < 1 (when H<sub>S</sub> is not 0, same issue with G<sub>ST</sub>). Nevertheless, to avoid a more serious situation in which heterozygosities of populations with one overlaped allele might be wrongly calculated if we only consider the two most frequent alleles, and at least one of the other alleles doesn't overlap between populations (e.g., P1: A=0.01, G=0.99; P2: A=0.02, C=0.98), we should still calculate heterozygosities through all existing allels (i.e. 1 - sum(square(p))).
4. In practical calculations, there might be sites where an allele is fixed or lost in both populations while not in other populations. This could lead to H<sub>T</sub> = 0, which violates the computation. Thus we will identify such cases and assign 0 to the pairwise F<sub>ST</sub>.
5. The start position of the first window of each contig is 1, and the following boundries of each window are midpoint between the current ending SNP and the next starting SNP. Beaware that this program is not strictly 'ungapped', because the most 3' SNP of a contig is not necessarily the end base of that contig, and also because some of the most 3' SNP windows are discarded due to the lack of accumulative heterozygosity.
6. The window is defined by accumulative heterozygosity (averaged among all samples analyzed) calculated from SNPs, which should have a summed heterozygosity >= 100 for each window and include ~10k sites.
7. The within-window F<sub>ST</sub> was calculated using weighted average across sites.
8. We get the across-window F<sub>ST</sub> as a weighted average by number of analyzed sites in each window.
9. The sample size used to calculate Fst is effective sample size, which is calculated at each site (instead of effective pool size, which is a component of calculating effective sample size, and is set as the actual pool size due to the lack of replicate population samples) based on [Ferretti et al, 2013](10.1111/mec.12522)
10. Keep in mind that F<sub>ST</sub> is a *relative* measure of population differentiation, and is strongly influenced by levels of within-population diversity.

## How we dealt with some tricky problems
1. There're many contigs within which the total number of SNPs is far less than what could accumulate enough heterozygosity to pass the threshold. On the other hand, we don't know the consecutive order among contigs. My solution is to skip these contigs for Fst calculation.
2. Similar to the above issue, even if a contig has enough SNPs to pass the heterozygosity threshold, there're likely remaining regions at the end of contigs that are not long enough to qualify for a 'window'. In this case, I skipped these tails for Fst calculation. Does this make sense?
3. A questions is whether to set the boundry of each window as the positions of starting and end SNPs? Here, we do not count filtered sites outside our windows (i.e. those before the first 5' window SNP and after the last 3' window SNP).