#!/bin/bash
# calculate the total length of unassigned contigs after correlation-based assignment
awk 'FILENAME==ARGV[1] && $0 !~ "#" {len[$7]=$10} \
FILENAME==ARGV[2] && $1 in len {tlen+=len[$1]} \
END{print tlen}' <(cat /Users/siyuansmac/test/GCA_013340165.1_LBDM_Dsuz_2.1.pri_assembly_report.txt) <(awk '$3 == "ambiguous" {print $1}' /Users/siyuansmac/bioinfo/project/suzukii_WGS/SV_analyses/CNV/detect_CNV/scripts/assignment_cor_amb.txt) 

# calculate the total length of unassigned contigs after assignment by Paris et al 2020
awk 'FILENAME==ARGV[1] && $0 !~ "#" {len[$7]=$10} \
FILENAME==ARGV[2] && $1 in len {tlen+=len[$1]} \
END{print tlen}' <(cat /Users/siyuansmac/test/GCA_013340165.1_LBDM_Dsuz_2.1.pri_assembly_report.txt) <(awk '$4 == "NA" {print $1}' /Users/siyuansmac/bioinfo/project/suzukii_WGS/SV_analyses/CNV/detect_CNV/scripts/assignment_cor_amb.txt) 
