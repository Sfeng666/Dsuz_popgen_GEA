#!/bin/bash
# calculate the mean depth of all contigs for correlation-based assignment

## 1. configure the conda enviroment
set -e
ENVDIR=WGS_analysis
. $ENVDIR/bin/activate

## 2. run the patch code
bam_list=bam_list.txt   # input table of bam files, with sample name at the begining of each '_' delimited string
genome_rep=*_assembly_report.txt    # input table of contig information from NCBI, containing mapping-based assignment of each contig 
mean_depth=mean_depth.txt   # output table of average depth of each contig in each sample

# ## 2.1. concatenate all Autosome-linked contigs and all X-linked contigs to caltulate the reference mean depth among all sites for each sample

# print the head line of sample names
awk -F "_" 'NR==1{printf "%s\t%s\t%s", "RefSeq-Accn", "Assigned-Molecule",$1} NR>1{printf "\t" $1} END{printf "\n"}' $bam_list > $mean_depth

# calculate the mean depth of all sites on Autosome-linked contigs
awk '$1 !~ /#/ && $3 != "na" && $3 != "X" {print $7}' $genome_rep \
|parallel zcat {}_mpileup.gz \
|cut -f4- \
|grep -E '[1-9]' \
|awk 'BEGIN{FS=OFS=" ";printf "%s\t%s", "auto-contigs", "autosome"} {for (i=1;i<=NF;i+=3) a[i]+=$i} END{for (i=1;i<=NF;i+=3) printf OFS a[i]/NR; printf "\n"}' >> $mean_depth

# calculate the mean depth of all sites on X-linked contigs
awk '$1 !~ /#/ && $3 == "X" {print $7}' $genome_rep \
|parallel zcat {}_mpileup.gz \
|cut -f4- \
|grep -E '[1-9]' \
|awk 'BEGIN{FS=OFS=" ";printf "%s\t%s", "X-contigs", "X"} {for (i=1;i<=NF;i+=3) a[i]+=$i} END{for (i=1;i<=NF;i+=3) printf OFS a[i]/NR; printf "\n"}' >> $mean_depth

## 2.2. calculate the mean depth of all sites on unplaced contigs (as well as placed contigs, to validate whether the correlation method is appropriate for assignment of contigs)
md() {
cont=$(echo $1|cut -d" " -f 1)
chr=$(echo $1|cut -d" " -f 2)
mean_depth=$2
awk -v cont=$cont -v chr=$chr 'BEGIN{FS=OFS=" ";printf "%s\t%s", cont, chr} {for (i=1;i<=NF;i+=3) a[i]+=$i} END{for (i=1;i<=NF;i+=3) printf OFS a[i]/NR; printf "\n"}' \
<(zcat $cont\_mpileup.gz |cut -f4- |grep -E '[1-9]') \
>> $mean_depth
}

export -f md
awk '$1 !~ /#/ {print $7, $3}' $genome_rep \
|parallel md {} $mean_depth
