#!/bin/bash

### 1. configure the conda enviroment ###
source /home/siyuan/biosoft/miniconda3/bin/activate wgs_suzukii

### 2. set parameters ###
rscript=/home/siyuan/jobs/suzukii_WGS/genetic_diff/ML_tree/plot/mincount_10_no_se_auto_X_over10/plot_ML_tree.R
wd=/home/siyuan/jobs/suzukii_WGS/genetic_diff/ML_tree/plot/mincount_10_no_se_auto_X_over10
dd=/home/siyuan/jobs/suzukii_WGS/genetic_diff/ML_tree/results/mincount_10_no_se_auto_X_over10
sp=/home/siyuan/jobs/suzukii_WGS/genetic_diff/ML_tree/treemix/treemix-1.13/src/plotting_funcs.R
cl=/home/siyuan/jobs/suzukii_WGS/genetic_diff/ML_tree/plot/pop_order_color_consistent.txt
log=run_plot_ML_tree.log

### 3. run the plotting R script ###
cd $wd
# Rscript $rscript $wd $dd $sp $cl winlen_500_mincov_12_mincount_1 >> $log 2>&1
Rscript $rscript $wd $dd $sp $cl winlen_500_mincov_12_mincount_10 >> $log 2>&1
# Rscript $rscript $wd $dd winlen_500_mincov_12_mincount_5 >> $log 2>&1