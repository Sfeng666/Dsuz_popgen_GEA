#!/bin/bash
# reformat the space-delimited file into a tab-delimited one
perl -nle '@a=split;print join "\t", @a' /Users/siyuansmac/bioinfo/temp/mean_depth.txt >/Users/siyuansmac/bioinfo/temp/mean_depth_format.txt

# generate a table with contig name, ref_seq accession and chromosomes assigned by Kapun
join -t "$(echo -e "\t")" <(cat /Users/siyuansmac/bioinfo/temp/GCF_013340165.1_LBDM_Dsuz_2.1.pri_assembly_report.txt | grep -v '#' | awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$7}' | sort -n -t"_" -k2) \
<(cat /Users/siyuansmac/bioinfo/temp/assignment_Kappun.txt | cut -f 1,2 | sort -n -t"_" -k2) | cut -f 2,3 > /Users/siyuansmac/bioinfo/temp/assignment_Kappun_refseq.txt
