#!/bin/bash
## this script is used for calculating nucleotide diversity for each sample (population), based on previously called SNPs, high-quality sites and indel positions.
## this is a modified version of 'calc_diversity.sh', with a change that we discard a whole site with at least one sample failing to meet the min- and max-cov filtering.

### script starts ###
echo start at time `date +%F'  '%H:%M`

### 1. configure the conda enviroment ###
source /opt/miniconda3/bin/activate WGS_analysis

### 2. set paths ###

## paths of parameters ##
wd=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/calc_diversity_maxcov_99/test_normalize/test_biallelic/rescale/auto_X_separated/calc_from_count/mincov_$1\_mincount_$2
log=$wd/calc_diversity.log
diversity=$wd/nucleotide_diversity.txt
temp=$wd/temp
sample_size=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/calc_diversity_maxcov_99/test_normalize/sample_size.txt

dir_snp=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/call_snp_maxcov_99/mincov_$1\_mincount_$2
snp_clean_contig=$dir_snp/{}_snp_clean.vcf
snp_clean=$dir_snp/snp_clean.vcf.gz
indel_positions_contig=$dir_snp/{1}_inDel-positions_20.txt

dir_mpileup=/home/siyuan/jobs/suzukii_WGS/06mpileup_rsync
mpileup_contig=$dir_mpileup/{1}_mpileup.gz
mincov=$1
maxcov=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/calc_maxcov_99/max_cov.txt
chr_assignment=w/assignment_cor_amb.txt
min_count=$2
min_freq=0
miss_frac=0.001
script=/home/siyuan/jobs/suzukii_WGS/calc_diversity_same_sites_ts/updated_poolsnp/filter_mpileup_by_cov.py

bam_list=$dir_mpileup/bam_list.txt

mkdir -p $temp

names=""
num_sample=0
while read -r line
do
    sample=$(basename $line _sort_merge_dedup_indel.bam)
    if [ $names ]
    then
        names+=",${sample}"
    else
        names+=$sample
    fi
    ((num_sample++))

done < $bam_list

### 3. run steps of the pipeline ###

echo -e "### Whole job starts ###\n" > $log

## define a function to calculate total heterozygosity for InDel-free SNPs within each contig for each sample

calc_hetero() {
local snp_clean_contig=$1
local log=$2
awk 'BEGIN{FS = "\t"} \
$0 !~ /#/ {i = 10; while (i <= NF) {split($i, attr, ":"); ref_ct = attr[2]; split(attr[3], ct, ","); for (j in ct) {ct_all[j] = ct[j]}; ct_all[ref] = ref_ct; \
max_ct=0; second_ct=0; for (ale in ct_all) {if (max_ct < ct_all[ale]) max_ct=ct_all[ale]}; \
for (ale in ct_all) {if (ct_all[ale] == max_ct) {delete ct_all[ale]; break}}; \
for (ale in ct_all) {if (second_ct < ct_all[ale]) second_ct=ct_all[ale]}; delete ct_all; \
total_ct = max_ct + second_ct; \
a[i] += (1 -  (max_ct/total_ct)^2 - (second_ct/total_ct)^2) * (total_ct)/(total_ct - 1); i++}} \
END{for (i = 10;i <= NF - 1;i++) printf "%f\t",a[i]; printf "%f\n",a[i]}' $snp_clean_contig

echo "heterozygosity calculated from $snp_clean_contig" >> $log
}

## Note: there're two ways to run shell command within awk and get the output back into awk, but the first seems to be ~20% faster than the second, using test data (/Users/siyuansmac/test/diversity_test/test.vcf). Nevertheless, it's much faster to use internal functions than invoking outside bash commands.
# awk 'BEGIN{FS=OFS=ORS="\t"; CMD="cut -d: -f5"} $0 !~ /#/ {i=10;while (i<=NF) {print $i |& CMD; close(CMD,"to"); CMD |& getline b; close(CMD); a[i]+=2*b*(1-b); i +=1}} END{for (i=10;i<=NF;i+=1) print a[i]; printf "\n"}' test.vcf
# awk 'BEGIN{FS=OFS=ORS="\t"} $0 !~ /#/ {i=10;while (i<=NF) {"echo " $i "|cut -d: -f5" | getline b; close("echo " $i "|cut -d: -f5"); a[i]+=2*b*(1-b); i+=1}} END{for (i=10;i<=NF;i+=1) print a[i]; printf "\n"}' test.vcf

## define a function to count the number of InDel-free and coverage-filtered sites
count_sites() {
local indel_positions_contig=$1
local mpileup_contig=$2
local min_cov=$3
local max_cov=$4
local min_count=$5
local min_freq=$6
local miss_frac=$7
local names=$8
local script=$9
local log=${10}
local temp=${11}


awk 'BEGIN{FS="\t"; sites=0} \
NR==FNR {D[$1$2]++; next} \
!($1$2 in D) && NR>FNR {print} ' \
<(cat $indel_positions_contig) \
<(zcat $mpileup_contig) \
| parallel \
--tmpdir $temp \
-k \
--pipe \
-j 100% \
--no-notice \
--cat python2.7 $script \
--mpileup  {} \
--min-cov $min_cov \
--max-cov $max_cov \
--min-count $min_count \
--min-freq $min_freq \
--miss-frac $miss_frac \
--names $names

echo "InDel-free and coverage-filtered sites counted from $mpileup_contig" >> $log
}

## For each contig,
## 3.1. calculate the accumulative value of heterozygosity for InDel-free SNPs in each sample, then output the values into a record file.
# note that the calculations are based on sites within autosomes-linked and X-linked (X-linked and ambiguous) contigs separately.

export -f calc_hetero

# autosome-linked contigs
echo -e "### Autosome-linked contigs: Calculation of accumulative heterozigosity starts ###\n" >> $log

{ time sum_heterozygosity_auto=$(awk '$4 == "auto-linked"{print $1}' <(join -1 1 -2 1 <(sort -k 1 $maxcov) <(sort -k 1 $chr_assignment)) \
| parallel --no-notice --tmpdir $temp -j 1000% \
calc_hetero $snp_clean_contig $log \
| awk 'BEGIN{FS="\t"} \
NR==2 {i=1; while (i<=NF) {samplesize[i]=$i; i++}} \
NR>3 {i=1;while (i<=NF) {a[i]+=$i; i++}} \
END{for (i=1;i<=NF-1;i+=1) printf "%f\t",a[i]*samplesize[i]/(samplesize[i] - 1); printf "%f\n",a[i]*samplesize[i]/(samplesize[i] - 1)}' <(cat $sample_size) - ) ;} 2>> $log

echo "Autosome-linked contigs: accumulative heterozygosity = $sum_heterozygosity_auto"  >> $log
echo -e "### Autosome-linked contigs: Calculation of accumulative heterozygosity ends ###\n" >> $log

# X-linked and ambiguous
echo -e "### X-linked contigs: Calculation of accumulative heterozigosity starts ###\n" >> $log

{ time sum_heterozygosity_X=$(awk '$4 == "X-linked" || $4 == "ambiguous"{print $1}' <(join -1 1 -2 1 <(sort -k 1 $maxcov) <(sort -k 1 $chr_assignment)) \
| parallel --no-notice --tmpdir $temp -j 1000% \
calc_hetero $snp_clean_contig $log \
| awk 'BEGIN{FS="\t"} \
NR==3 {i=1; while (i<=NF) {samplesize[i]=$i; i++}} \
NR>3 {i=1;while (i<=NF) {a[i]+=$i; i++}} \
END{for (i=1;i<=NF-1;i+=1) printf "%f\t",a[i]*samplesize[i]/(samplesize[i] - 1); printf "%f\n",a[i]*samplesize[i]/(samplesize[i] - 1)}' <(cat $sample_size) - ) ;} 2>> $log

echo "X-linked contigs: accumulative heterozygosity = $sum_heterozygosity_X"  >> $log
echo -e "### X-linked contigs: Calculation of accumulative heterozygosity ends ###\n" >> $log

## 3.2. count the number of InDel-free filtered SNPs, then output the values into a record file

# autosome-linked contigs
echo -e "### Autosome-linked contigs: Counting of InDel-free filtered SNPs starts ###\n" >> $log

{ time num_snps_auto=$(awk '$4 == "auto-linked"{print $1}' <(join -1 1 -2 1 <(sort -k 1 $maxcov) <(sort -k 1 $chr_assignment)) \
| parallel --no-notice -j 1000% \
zgrep -v "^#" $snp_clean_contig \
| wc -l \
| awk -v num_sample=$num_sample '{sites=$1} \
END{for (i=1;i<=num_sample - 1;i+=1) printf "%d\t",sites; printf "%d\n",sites}') ;} 2>> $log

echo "Autosome-linked contigs: number of InDel-free filtered SNPs = $num_snps_auto"  >> $log
echo -e "### Autosome-linked contigs: Counting of InDel-free filtered SNPs ends ###\n" >> $log

# X-linked and ambiguous
echo -e "### X-linked contigs: Counting of InDel-free filtered SNPs starts ###\n" >> $log

{ time num_snps_X=$(awk '$4 == "X-linked" || $4 == "ambiguous"{print $1}' <(join -1 1 -2 1 <(sort -k 1 $maxcov) <(sort -k 1 $chr_assignment)) \
| parallel --no-notice -j 1000% \
zgrep -v "^#" $snp_clean_contig \
| wc -l \
| awk -v num_sample=$num_sample '{sites=$1} \
END{for (i=1;i<=num_sample - 1;i+=1) printf "%d\t",sites; printf "%d\n",sites}') ;} 2>> $log

echo "X-linked contigs: number of InDel-free filtered SNPs = $num_snps_X"  >> $log
echo -e "### X-linked contigs: Counting of InDel-free filtered SNPs ends ###\n" >> $log


## 3.3. count the number of InDel-free filtered sites, then output the values into a record file

# autosome-linked contigs
echo -e "### Autosome-linked contigs: Counting of InDel-free filtered sites starts ###\n" >> $log

export -f count_sites
{ time num_sites_auto=$(awk 'BEGIN{OFS="\t"} $4 == "auto-linked"{print $1, $2}' <(join -1 1 -2 1 <(sort -k 1 $maxcov) <(sort -k 1 $chr_assignment)) \
| parallel --no-notice --tmpdir $temp -j 1000% --colsep '\t' \
count_sites $indel_positions_contig $mpileup_contig $mincov $maxcov $min_count $min_freq $miss_frac $names $script $log $temp \
| awk -v num_sample=$num_sample 'BEGIN{sites=0} \
{sites+=$1} \
END{for (i=1;i<=num_sample - 1;i+=1) printf "%d\t",sites; printf "%d\n",sites}') ;} 2>> $log

echo "Autosome-linked contigs: number of InDel-free filtered sites = $num_sites_auto"  >> $log
echo -e "### Autosome-linked contigs: Counting of InDel-free filtered sites ends ###\n" >> $log

# X-linked and ambiguous
echo -e "### X-linked contigs: Counting of InDel-free filtered sites starts ###\n" >> $log

export -f count_sites
{ time num_sites_X=$(awk 'BEGIN{OFS="\t"} $4 == "X-linked" || $4 == "ambiguous"{print $1, $2}' <(join -1 1 -2 1 <(sort -k 1 $maxcov) <(sort -k 1 $chr_assignment)) \
| parallel --no-notice --tmpdir $temp -j 1000% --colsep '\t' \
count_sites $indel_positions_contig $mpileup_contig $mincov $maxcov $min_count $min_freq $miss_frac $names $script $log $temp \
| awk -v num_sample=$num_sample 'BEGIN{sites=0} \
{sites+=$1} \
END{for (i=1;i<=num_sample - 1;i+=1) printf "%d\t",sites; printf "%d\n",sites}') ;} 2>> $log

echo "X-linked contigs: number of InDel-free filtered sites = $num_sites_X"  >> $log
echo -e "### X-linked contigs: Counting of InDel-free filtered sites ends ###\n" >> $log

## 3.3. Estimate nucleotide diversity as average heterozygosity, with the accumulative heterozygosity as numerator, and with the number of InDel-free filtered sites as denominator

awk 'BEGIN{FS="_"} NR==1{printf $1} NR>1{printf "\t%s", $1} END{printf "\n"}' $bam_list > $diversity

# autosome-linked contigs
echo -e "### Autosome-linked contigs: Estimation of nucleotide diversity starts ###\n" >> $log

awk 'NR==1{i=1; while (i<=NF) {sh[i]=$i; if (i==1) printf $i; else printf "\t%f", $i; i++}; printf "\n"} \
NR==2{i=1; while (i<=NF) {nsnp[i]=$i; if (i==1) printf $i; else printf "\t%d", $i; i++}; printf "\n"}
NR==3{i=1; while (i<=NF) {ns[i]=$i; if (i==1) printf $i; else printf "\t%d", $i; i++}; printf "\n"} \
END{i=1; while (i<=NF) {if (i==1) printf sh[i]/ns[i]; else printf "\t%f", sh[i]/ns[i]; i++}; printf "\n"}' \
<(echo $sum_heterozygosity_auto) <(echo $num_snps_auto) <(echo $num_sites_auto) >> $diversity

echo -e "### Autosome-linked contigs: Estimation of nucleotide diversity ends ###\n" >> $log

# X-linked and ambiguous
echo -e "### X-linked contigs: Estimation of nucleotide diversity starts ###\n" >> $log

awk 'NR==1{i=1; while (i<=NF) {sh[i]=$i; if (i==1) printf $i; else printf "\t%f", $i; i++}; printf "\n"} \
NR==2{i=1; while (i<=NF) {nsnp[i]=$i; if (i==1) printf $i; else printf "\t%d", $i; i++}; printf "\n"}
NR==3{i=1; while (i<=NF) {ns[i]=$i; if (i==1) printf $i; else printf "\t%d", $i; i++}; printf "\n"} \
END{i=1; while (i<=NF) {if (i==1) printf sh[i]/ns[i]; else printf "\t%f", sh[i]/ns[i]; i++}; printf "\n"}' \
<(echo $sum_heterozygosity_X) <(echo $num_snps_X) <(echo $num_sites_X) >> $diversity

echo -e "### X-linked contigs: Estimation of nucleotide diversity ends ###\n" >> $log

echo -e "### Whole job done ###\n" >> $log

rm -r $temp

### script finishes ###
echo finish at time `date +%F'  '%H:%M`