#! /bin/bash
# run a modified GSEA script to generate list of outlier-associated genes
cd /home/siyuan/jobs/suzukii_WGS/EAA/downstream/GSEA/random_SNP/scripts/output_gene_list/all_trimmed_outliers/scripts

for env_var in Altitude Annual_mean_wind_speed Annual_Precipitation Precipitation_Seasonality Mean_Diurnal_Range Mean_Temperature_of_Coldest_Quarter Mean_Temperature_of_Warmest_Quarter Ratio_Built_to_Vegetation Ratio_Crop_to_Forest
do
    bash run_GSEA.sh $env_var &
done