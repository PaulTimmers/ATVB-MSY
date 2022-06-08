#!/bin/bash

#----
# QC variants
#----

MAIN=/exports/igmm/eddie/wilson-lab/projects/prj_180_ygen_ukbb
mkdir -p ${MAIN}/p001_qc_variants
cd ${MAIN}/p001_qc_variants

# Variants within ISOGG
awk 'NR == FNR {pos[$3]++; next} pos[$4] > 0 {print $2}' \
../data/d001_isogg_snps_hg19.tsv  ../data/d002_ukbb_snps.hg19.bim \
> st001_01_shared_variants.txt


# Perform variant quality control
sample_qc=../data/d003_ukbb_sample_qc.tsv
min_snp_call_rate=0.95

plink2 \
--bfile ../data/d002_ukbb_snps.hg19 \
--extract st001_01_shared_variants.txt \
--keep <(awk -v infer_sex_col=12 '$infer_sex_col == "M" {print $1,$1}' ${sample_qc}) \
--geno `echo 1-$min_snp_call_rate | bc` \
--mac 1 \
--make-bed \
--out st001_02_all_samples
 

#----
# QC samples
#----

sample_qc=../data/d003_ukbb_sample_qc.tsv
sample_rel=../data/d003_ukbb_rel.dat

# Males
awk -v infer_sex_col=12 \
'BEGIN{OFS="\t"} $infer_sex_col == "M" {n++; print $1 > "st001_03_iids.txt"} END{print "Males",n}' ${sample_qc} > st001_03_stats.txt

# Heterozygosity and missingess outliers
awk -v het_mis_col=20 \
'BEGIN{OFS="\t"} NR == FNR {iid[$1]++; next} iid[$1] > 0 && $het_mis_col == 0 {n++; print $1 > "st001_03_iids.txt"} END{print "Outliers in heterozygosity and missingness",n}' \
st001_03_iids.txt ${sample_qc} >> st001_03_stats.txt

# Affymetrix QC metrics
awk -v cluster_col=8 -v cluster_thresh=97 -v dqc_col=9 -v dqc_thresh=0.82 \
'BEGIN{OFS="\t"} NR == FNR {iid[$1]++; next} iid[$1] > 0 && $cluster_col >= cluster_thresh && $dqc_col >= dqc_thresh {n++; print $1 > "st001_03_iids.txt"} END{print "Affymetrix QC metrics",n}' \
st001_03_iids.txt ${sample_qc} >> st001_03_stats.txt

# Missingness check
awk -v missing_col=17 -v missing_thresh=0.02 \
'BEGIN{OFS="\t"} NR == FNR {iid[$1]++; next} iid[$1] > 0 && $missing_col <= missing_thresh {n++; print $1 > "st001_03_iids.txt"} END{print "Missingness check",n}' \
st001_03_iids.txt ${sample_qc} >> st001_03_stats.txt

# Sex inconsistency
awk  -v sex_aneu_col=21 \
'BEGIN{OFS="\t"} NR == FNR {iid[$1]++; next} iid[$1] > 0 && $sex_aneu_col == 0  {n++; print $1 > "st001_03_iids.txt"} END{print "Sex inconsistency",n}' \
st001_03_iids.txt ${sample_qc} >> st001_03_stats.txt


# Relatedness
Rscript ../scripts/a001_03_remove_relatives.R


# White British ancestry
awk -v white_briti_col=25 \
'BEGIN{OFS="\t"} NR == FNR {iid[$1]++; next} iid[$1] > 0 && $white_briti_col == 1  {n++; print $1 > "st001_03_iids.txt"} END{print "White British ancestry",n}' \
st001_03_iids.txt ${sample_qc} >> st001_03_stats.txt

