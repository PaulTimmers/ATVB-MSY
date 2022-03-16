#!/bin/bash

#----
# Setup environment
#----

MAIN=/exports/igmm/eddie/wilson-lab/projects/prj_180_ygen_ukbb
cd ${MAIN}


#----
# Create haplogroups
#----


mkdir -p ${MAIN}/p002_yhaplo
cd ${MAIN}/p002_yhaplo


# Input
#------

mkdir -p ${MAIN}/p002_yhaplo/converted

# Create temporary MAP & PED files
plink \
--bfile ../p001_qc_variants/st001_02_all_samples \
--recode \
--out converted/st002_01_all_samples && rm ../p001_qc_variants/st001_02_all_samples.*


# Convert PED input to GENOS
source ../scripts/venv/bin/activate
convert_to_genos converted/st002_01_all_samples.ped && rm converted/st002_01_all_samples.ped converted/st002_01_all_samples.map
mv converted/ input/


# Run yhaplo
qsub -N yhaplo -l h_vmem=32G -l h_rt=48:00:00 -j y -o log_yhaplo.log -cwd -V <<"QSUB"
truncate -s 0 log_yhaplo.log
source ../scripts/venv/bin/activate
yhaplo -i input/st002_01_all_samples.genos.txt
QSUB


# Convert to wilson-format
R --quiet --no-save <<"Rcode"
df <- data.table::fread("output/haplogroups.st002_01_all_samples.txt", header=FALSE, col.names=c("iid.iid","h1","haplogroup_short","haplogroup"))
df <- df[,.(iid=gsub("-.*$","",iid.iid), haplogroup, haplogroup_short)][order(iid),]
df <- df[iid > 0,]
data.table::fwrite(df, "st002_02_haplogroups.tsv", sep="\t", na="NA", quote=FALSE)
Rcode

# Count
echo -e "haplogroup\tn\tpct" > st002_02_haplogroups_count.tsv
awk -v OFS="\t" 'NR > 1 {sum[$3]++; total++} END{for(haplogroup in sum){print haplogroup, sum[haplogroup], sum[haplogroup]/total}}' \
st002_02_haplogroups.tsv | sort -k2gr,2 >> st002_02_haplogroups_count.tsv



#----
# Sanity check
#----

# Get ethnicities
R --quiet --no-save <<Rcode
ph <- data.table::fread("../data/d003_ukbb_data.tsv", select=c("f.eid","f.21000.0.0"), col.names=c("iid","ethnicity"))
data.table::fwrite(ph, "st002_03_ethnicity.tsv", sep="\t", quote=FALSE, na="NA")
Rcode

# List most common haplogroup by major ethnic groups
Rscript ../scripts/a002_04_most_common_by_ethnicity.R
