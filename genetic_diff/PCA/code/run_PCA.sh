#!/bin/bash

### 1. configure the conda enviroment ###
source /opt/miniconda3/bin/activate WGS_analysis

### 2. set parameters ###
rscript=/home/siyuan/jobs/suzukii_WGS/genetic_diff/PCA/plot/plot_PCA.R
wd=/home/siyuan/jobs/suzukii_WGS/genetic_diff/PCA/plot
dd=/home/siyuan/jobs/suzukii_WGS/genetic_diff/PCA/results

### 3. run the plotting R script ###
Rscript $rscript $wd $dd mincov_12_mincount_1 >> run_PCA_maf.log 2>&1
# Rscript $rscript $wd $dd mincov_12_mincount_10 >> run_PCA_maf.log 2>&1
# Rscript $rscript $wd $dd mincov_12_mincount_20 >> run_PCA_maf.log 2>&1
# Rscript $rscript $wd $dd mincov_12_mincount_3 >> run_PCA_maf.log 2>&1
# Rscript $rscript $wd $dd mincov_12_mincount_5 >> run_PCA_maf.log 2>&1
# Rscript $rscript $wd $dd mincov_12_mincount_50 >> run_PCA_maf.log 2>&1
# Rscript $rscript $wd $dd mincov_12_mincount_8 >> run_PCA_maf.log 2>&1
# Rscript $rscript $wd $dd mincov_20_mincount_10 >> run_PCA_maf.log 2>&1