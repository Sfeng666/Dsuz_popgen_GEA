#! /bin/bash

cd /home/sfeng77/jobs/EAA/bayescenv/chtc/split_jobs/test_para_p

## this script submit dag for all environmental variables
# for env_var in Altitude 
for env_var in Altitude  Annual_Precipitation Precipitation_Seasonality Mean_Diurnal_Range Mean_Temperature_of_Coldest_Quarter Mean_Temperature_of_Warmest_Quarter Ratio_Built_to_Vegetation
# for env_var in Altitude Annual_mean_wind_speed Annual_Precipitation Precipitation_Seasonality Mean_Diurnal_Range Mean_Temperature_of_Coldest_Quarter Mean_Temperature_of_Warmest_Quarter Ratio_Built_to_Vegetation
do
    for p in 0.9
    do 
        condor_submit_dag results/$env_var\_$p/split_jobs_$env_var\_$p.dag
    done
done