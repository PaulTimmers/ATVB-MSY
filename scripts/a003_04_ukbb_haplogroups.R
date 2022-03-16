#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------------------------------------
# Setup environment
#----

source("../scripts/print_script_name.R")
options(width=200)
options(digits=6)

# R packages

library(data.table)
library(readxl)

# Variables
min_cases <- 40
min_sample <- 100


# Functions
heading <- function(sentence) {
	writeLines(paste("\n\n\n=======================\n\n",
	sentence,"\n==========================="))
}


#----------------------------------------------------------------------------------------------------------------
# Start
#----

sessionInfo()


#----
# Load data
#----

ph <- fread("st003_01_ukbb_base_data.tsv.gz")
biochem <- fread("st003_03_ukbb_biochem_data.tsv.gz")
haplogroups <- fread("../p002_yhaplo/st002_02_haplogroups.tsv")

# Merge
ph <- biochem[ph,,on="iid"]
ph <- haplogroups[ph,,on="iid"]


#----
# Haplogroups
#----

heading("y haplogroups")

writeLines(paste0("\nIndividuals without haplogroup: ",nrow(ph[is.na(ph$haplogroup) | ph$haplogroup == "A", ])))
ph <- ph[!(is.na(ph$haplogroup) | ph$haplogroup == "A"), ]



haplomap <- fread("../p002_yhaplo/output/haplogroups.st002_01_all_samples.txt", header=FALSE, col.names=c("iid","del","short","long"))
haplomap <- haplomap[!duplicated(paste0(short,long)),.(long,short)]

counts <- ph[,.(N_total=.N,
	N_cvd_cases=sum(cvd,na.rm=T),
	N_cvd_controls=sum(!cvd, na.rm=T), 
	N_sys=sum(!is.na(sys_blood_pressure)),
	N_dia=sum(!is.na(dia_blood_pressure)), 
	N_cholesterol=sum(!is.na(total_cholesterol)),
	N_ldl=sum(!is.na(ldl_cholesterol)),
	N_hdl=sum(!is.na(hdl_cholesterol)),
	N_triglyc=sum(!is.na(log_triglycerides))),by="haplogroup"]

counts <- haplomap[counts,,on=c(long="haplogroup")]
counts <- counts[order(long)]
fwrite(counts, "st003_04_pheno_counts.csv", sep=",", na="NA", quote=FALSE)


#----
# Groupings
#----

# Haplogroup names were fixed and new, additional (sub)groups were created based on a minimum of 40 cases
# These data were saved in st003_04_haplo_groupings.xlsx

# Load fixed names
ht_df <- data.table(read_xlsx("st003_04_haplo_groups.xlsx", sheet="Haplogroups", skip=1))
names(ht_df)[c(1,2,4,5)] <- c("haplogroup","haplogroup_short","wilson_haplo","wilson_haplo_short")
ht_df <- ht_df[!haplogroup=="Total"]

# Fix names
ph$haplogroup <- ht_df[match(ph$haplogroup_short, haplogroup_short), wilson_haplo]
ph$haplogroup_short <- ht_df[match(ph$haplogroup_short, haplogroup_short),wilson_haplo_short]


#----
# Export
#----

writeLines("\nCounts:")
summary <- sort(with(ph, table(haplogroup_short)), decreasing=TRUE)
summary[summary < 20] <- "<20"
print(summary)

variables <- c("iid", "haplogroup", "haplogroup_short")
ph1 <- ph[,.SD,.SDcols=variables]
fwrite(ph1, "st003_04_ukbb_haplogroups.tsv.gz", sep="\t", quote=TRUE, na="NA", compress="gzip")