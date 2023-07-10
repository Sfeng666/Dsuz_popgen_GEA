#!/bin/bash

### 1. configure the conda enviroment ###
source /opt/miniconda3/bin/activate WGS_analysis

### 2. set parameters ###
rscript=/home/siyuan/jobs/suzukii_WGS/genetic_diff/window_diversity/plot/plot_win_diversity.R
wd=/home/siyuan/jobs/suzukii_WGS/genetic_diff/window_diversity/plot
dd=/home/siyuan/jobs/suzukii_WGS/genetic_diff/window_diversity/results

### 3. run the plotting R script ###
Rscript $rscript $wd $dd winlen_125000_mincov_12_mincount_1 >> run_plot_win_diversity.log 2>&1