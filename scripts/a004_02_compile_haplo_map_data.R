#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------------------------------------
# Setup environment
#----

options(width=200)
options(digits=6)

# Data management
library(data.table)
library(readxl)
library(dplyr)

# Plotting
library(ggplot2)
library(viridis)

# Static
min_cases <- 40
min_sample <- 100
min_n_per_region <- 100

sessionInfo()

#----------------------------------------------------------------------------------------------------------------
# START
#----

ht_df <- data.table(read_xlsx("../p003_phenotypes/st003_04_haplo_groups.xlsx", sheet="Haplogroups", skip=1))
names(ht_df)[c(1,2,4,5)] <- c("haplogroup","haplogroup_short","wilson_haplo","wilson_haplo_short")
ht_df <- ht_df[!haplogroup=="Total"]

# Load groupings
htg_df <- data.table(read_xlsx("../p003_phenotypes/st003_04_haplo_groups.xlsx", sheet="Grouping"))
names(htg_df)[c(1:2)] <- c("wilson_haplo","wilson_haplo_short")
htg_df <- htg_df[!wilson_haplo=="Total"]

# Load group members
htgrouping <- data.table(read_xlsx("../p003_phenotypes/st003_04_haplo_groups.xlsx", sheet="Group Members"))



#----
# Spatial coordinate data
#----

map <- fread("st004_01_coordinate_data.tsv.gz")

# Load hierarchy
region_mapper <- fread("../data/d003_ukbb_geography/all-areas-lookup.csv")

# Format hierarchy mapper
region_mapper1 <- region_mapper[,.(place_birth=ifelse(CTRY_NAME=="England",RGN_GSS_ID,CTRY_GSS_ID),place_birth2=LA_GSS_ID,place_birth3=WED_GSS_ID)][!duplicated(paste0(place_birth,place_birth2,place_birth3))]

# Format country mapper
country_name_mapper <- region_mapper[,.(place_birth=CTRY_NAME,id=CTRY_GSS_ID)]
country_name_mapper <- country_name_mapper[!duplicated(paste0(place_birth, id))]


#----
# Load UK Biobank data
#----

# Load file matching iids with regions
ukbb_regions <- fread("../p003_phenotypes/st003_02_ukbb_geo_data.tsv.gz", select=c("iid","place_birth","rgn_birth", "dist_birth", "ward_birth"))

# Rename place birth to ids in geography data
ukbb_regions$place_birth <- country_name_mapper[match(ukbb_regions$place_birth, place_birth), id]

# Create variables matching region names
ukbb_regions[,place_birth:=ifelse(is.na(rgn_birth), place_birth, rgn_birth)]
ukbb_regions[,place_birth2:=ifelse(is.na(dist_birth), place_birth, dist_birth)]
ukbb_regions[,place_birth3:=ifelse(is.na(ward_birth), place_birth2, ward_birth)]

# Load haplogroups
ph <- fread("../p002_yhaplo/st002_02_haplogroups.tsv", data.table=FALSE)
ph$haplogroup <- ht_df[match(ph$haplogroup_short, haplogroup_short), wilson_haplo]
ph$haplogroup_short <- ht_df[match(ph$haplogroup_short, haplogroup_short),wilson_haplo_short]

# Load list of unrelated, white British individuals
iid_df <- fread("../p001_qc_variants/st001_03_iids.txt", col.names="iid")

# Exclude non-white, related individuals from main dataframe
ph[!ph$iid %in% iid_df$iid, "haplogroup"] <- NA
ph[!ph$iid %in% iid_df$iid, "haplogroup_short"] <- NA


#----
# Create MSY vars
#----

haplos <- unique(c(ht_df$wilson_haplo_short, htgrouping$Group))

# Define groups based on minimum cases
cvd_groups <- c(ht_df[Cases >= min_cases, wilson_haplo_short], 
  htg_df[`CVD Cases` >= min_cases, wilson_haplo_short])

# Define groups based on minimum sample
quant_groups <- c(ht_df[Count >= min_sample, wilson_haplo_short], 
  htg_df[Count >= min_sample, wilson_haplo_short])



#----
# Create fill variable
#----

# Merge MSY data with region data
df_phenotypes1 <- merge(ukbb_regions, ph, on="iid")
df_phenotypes1 <- df_phenotypes1[,.(iid, haplogroup_short, place_birth, place_birth2, place_birth3)]

# Create additional entries for grouped haplogroups
df_phenotypes2 <- copy(df_phenotypes1)
for(haplo in unique(htgrouping$Group)) {
  dat <- df_phenotypes1[haplogroup_short %in% htgrouping[Group==haplo, Members], ]
  dat[,haplogroup_short:=haplo]
  df_phenotypes2 <- rbind(df_phenotypes2, dat)
}


haplos <- unique(c(cvd_groups, quant_groups))

# Summarise
#----------

# Total per region

sample_df1 <- df_phenotypes2[!duplicated(iid),.(TOT1=.N),by=c("place_birth")]
sample_df2 <- df_phenotypes2[!duplicated(iid),.(TOT2=.N),by=c("place_birth2")]
sample_df3 <- df_phenotypes2[!duplicated(iid),.(TOT3=.N),by=c("place_birth3")]

regions_sample1 <- merge(region_mapper1,sample_df1,by="place_birth", all=T)[!is.na(place_birth3),.(place_birth3, TOT1)]
regions_sample2 <- merge(region_mapper1,sample_df2,by="place_birth2", all=T)[!is.na(place_birth3),.(place_birth3, TOT2)]
regions_sample3 <- merge(region_mapper1,sample_df3,by="place_birth3", all=T)[!is.na(place_birth3),.(place_birth3, TOT3)]

regions_sample_df <- merge(regions_sample3,regions_sample2,by=c("place_birth3"), all=TRUE) 
regions_sample_df <- merge(regions_sample_df,regions_sample1,by=c("place_birth3"), all=TRUE) 
regions_sample_df <- regions_sample_df[!duplicated(paste0(place_birth3, TOT1, TOT2, TOT3)),]

regions_sample_df <- data.table(regions_sample_df, haplogroup_short=rep(haplos,length(unique(regions_sample_df$place_birth3))))


# N haplos per region

summary_df1 <- df_phenotypes2[,.(N1=.N),by=c("haplogroup_short", "place_birth")][!is.na(haplogroup_short)]
summary_df2 <- df_phenotypes2[,.(N2=.N),by=c("haplogroup_short", "place_birth2")][!is.na(haplogroup_short)]
summary_df3 <- df_phenotypes2[,.(N3=.N),by=c("haplogroup_short", "place_birth3")][!is.na(haplogroup_short)]

regions_summary1 <- merge(region_mapper1,summary_df1,by="place_birth",allow.cartesian=TRUE, all=TRUE)[!is.na(place_birth3),.(place_birth3, haplogroup_short, N1)]
regions_summary2 <- merge(region_mapper1,summary_df2,by="place_birth2", allow.cartesian=TRUE, all=TRUE)[!is.na(place_birth3),.(place_birth3, haplogroup_short, N2)]
regions_summary3 <- merge(region_mapper1,summary_df3,by="place_birth3", allow.cartesian=TRUE, all=TRUE)[!is.na(place_birth3),.(place_birth3, haplogroup_short, N3)]

regions_summary_df <- merge(regions_summary3,regions_summary2,by=c("place_birth3","haplogroup_short"), all=TRUE) 
regions_summary_df <- merge(regions_summary_df,regions_summary1,by=c("place_birth3","haplogroup_short"), all=TRUE) 
regions_summary_df <- regions_summary_df[!duplicated(paste0(place_birth3,haplogroup_short, N1, N2, N3)),]


# Combine n_haplogroup and n_total data
summary_sample_df <- merge(regions_summary_df, regions_sample_df, by=c("place_birth3", "haplogroup_short"), allow.cartesian=TRUE, all=T)
summary_sample_df[is.na(N1), N1:=0]


# Fill missing values with higher-order region values
summary_sample_df[TOT1 >= min_n_per_region & !is.na(N1) & !is.na(TOT1), PCT:=N1/TOT1]
summary_sample_df[TOT2 >= min_n_per_region & !is.na(N2) & !is.na(TOT2), PCT:=N2/TOT2]
summary_sample_df[TOT3 >= min_n_per_region & !is.na(N3) & !is.na(TOT3), PCT:=N3/TOT3]








#----
# Export data for plotting
#----

# Format for plotting
plot_df <- summary_sample_df[!is.na(haplogroup_short) & place_birth3 %in% map$id, .(id=place_birth3, haplogroup_short, pct_haplo=PCT)]
plot_df <- plot_df[!duplicated(paste0(id, haplogroup_short)),]


# Select haplos to plot
common_haplos <- sort(unique(c(cvd_groups, quant_groups)))


# Export data
fwrite(plot_df, "st004_02_map_plot_data.tsv.gz", compress="gzip", sep="\t", na="NA")
fwrite(data.table(common_haplos), "st004_02_map_plot_haplos.txt", sep="\t", na="NA", quote=FALSE, col.names=FALSE)
