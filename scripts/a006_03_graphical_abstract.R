#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------------------------------------
# Setup environment
#----

options(width=200)
options(digits=6)
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(viridis)

# Static
min_n_per_region <- 100


#----
# Load data
#----

map <- fread("../p004_maps/st004_01_coordinate_data.tsv.gz")
plot_df <- fread("../p004_maps/st004_02_map_plot_data.tsv.gz")[haplogroup_short=="I1-M253",.(id,pct_haplo)]
plot_df[,pct_haplo_bin:=as.numeric(cut(pct_haplo,breaks=unique(c(-Inf,quantile(pct_haplo, probs=1:19/20,na.rm=T),Inf))))]



#----
# Summarise CAD by region
#----



# Load hierarchy
region_mapper <- fread("../data/d003_ukbb_geography/all-areas-lookup.csv")

# Format hierarchy mapper
region_mapper1 <- region_mapper[,.(place_birth=ifelse(CTRY_NAME=="England",RGN_GSS_ID,CTRY_GSS_ID),place_birth2=LA_GSS_ID,place_birth3=WED_GSS_ID)][!duplicated(paste0(place_birth,place_birth2,place_birth3))]

# Format country mapper
country_name_mapper <- region_mapper[,.(place_birth=CTRY_NAME,id=CTRY_GSS_ID)]
country_name_mapper <- country_name_mapper[!duplicated(paste0(place_birth, id))]


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

# Load CAD phenotypes
ph  <- fread("../p003_phenotypes/st003_01_ukbb_base_data.tsv.gz", select=c("iid","cvd"))

# Load list of unrelated, white British individuals
iid_df <- fread("../p001_qc_variants/st001_03_iids.txt", col.names="iid")

# Exclude non-white, related individuals from main dataframe
ph[!ph$iid %in% iid_df$iid, "cvd"] <- NA


# Create fill variable
#----

# Merge MSY data with region data
df_phenotypes1 <- merge(ukbb_regions, ph, on="iid")
df_phenotypes1 <- df_phenotypes1[,.(iid, cvd, place_birth, place_birth2, place_birth3)]


# Summarise
#----------

# Total per region

sample_df1 <- df_phenotypes1[!duplicated(iid) & !is.na(cvd),.(TOT1=.N),by=c("place_birth")]
sample_df2 <- df_phenotypes1[!duplicated(iid) & !is.na(cvd),.(TOT2=.N),by=c("place_birth2")]
sample_df3 <- df_phenotypes1[!duplicated(iid) & !is.na(cvd),.(TOT3=.N),by=c("place_birth3")]

regions_sample1 <- merge(region_mapper1,sample_df1,by="place_birth", all=T)[!is.na(place_birth3),.(place_birth3, TOT1)]
regions_sample2 <- merge(region_mapper1,sample_df2,by="place_birth2", all=T)[!is.na(place_birth3),.(place_birth3, TOT2)]
regions_sample3 <- merge(region_mapper1,sample_df3,by="place_birth3", all=T)[!is.na(place_birth3),.(place_birth3, TOT3)]

regions_sample_df <- merge(regions_sample3,regions_sample2,by=c("place_birth3"), all=TRUE) 
regions_sample_df <- merge(regions_sample_df,regions_sample1,by=c("place_birth3"), all=TRUE) 
regions_sample_df <- regions_sample_df[!duplicated(paste0(place_birth3, TOT1, TOT2, TOT3)),]


# N cvd per region

summary_df1 <- df_phenotypes1[!is.na(cvd),.(N1=sum(cvd)),by=c("place_birth")]
summary_df2 <- df_phenotypes1[!is.na(cvd),.(N2=sum(cvd)),by=c("place_birth2")]
summary_df3 <- df_phenotypes1[!is.na(cvd),.(N3=sum(cvd)),by=c("place_birth3")]

regions_summary1 <- merge(region_mapper1,summary_df1,by="place_birth",allow.cartesian=TRUE, all=TRUE)[!is.na(place_birth3),.(place_birth3, N1)]
regions_summary2 <- merge(region_mapper1,summary_df2,by="place_birth2", allow.cartesian=TRUE, all=TRUE)[!is.na(place_birth3),.(place_birth3, N2)]
regions_summary3 <- merge(region_mapper1,summary_df3,by="place_birth3", allow.cartesian=TRUE, all=TRUE)[!is.na(place_birth3),.(place_birth3, N3)]

regions_summary_df <- merge(regions_summary3,regions_summary2,by=c("place_birth3"), all=TRUE) 
regions_summary_df <- merge(regions_summary_df,regions_summary1,by=c("place_birth3"), all=TRUE) 
regions_summary_df <- regions_summary_df[!duplicated(paste0(place_birth3, N1, N2, N3)),]


# Combine n_haplogroup and n_total data
summary_sample_df <- merge(regions_summary_df, regions_sample_df, by=c("place_birth3"), allow.cartesian=TRUE, all=T)
summary_sample_df[is.na(N1), N1:=0]


# Fill missing values with higher-order region values
summary_sample_df[TOT1 >= min_n_per_region & !is.na(N1) & !is.na(TOT1), PCT:=N1/TOT1]
summary_sample_df[TOT2 >= min_n_per_region & !is.na(N2) & !is.na(TOT2), PCT:=N2/TOT2]
summary_sample_df[TOT3 >= min_n_per_region & !is.na(N3) & !is.na(TOT3), PCT:=N3/TOT3]

# Format for plotting
cad_plot_df <- summary_sample_df[place_birth3 %in% map$id, .(id=place_birth3, pct_cad=PCT)]
cad_plot_df <- cad_plot_df[!duplicated(id),]
cad_plot_df[,pct_cad_bin:=as.numeric(cut(pct_cad,breaks=unique(c(-Inf,quantile(pct_cad, probs=1:19/20),Inf))))]


#----
# Create individual plots
#----

sp_map_dat <- plot_df[map,,on="id"]
sp_map_cad <- cad_plot_df[map,,on="id"]
	
p1 <- ggplot(sp_map_dat) +
geom_polygon(aes(x = long, y = lat, group = group, fill = pct_haplo_bin), color=NA, size=0.1) +
scale_fill_viridis(guide="none") +
coord_fixed(xlim=c(75000, 750000)) + # Reduce horizontal whitespace around UK
theme_void() + 
theme(legend.title=element_text(size=14, vjust=1), legend.text=element_text(size=12))
ggsave("st006_03_graphical_abstract_i1.png", p1, width=7.5, height=7.5, dpi=300, units="cm")


p2 <- ggplot(sp_map_cad) +
geom_polygon(aes(x = long, y = lat, group = group, fill = pct_cad_bin), color=NA, size=0.1) +
scale_fill_gradient(low="#fff5f0", high="#67000d", guide="none") +
coord_fixed(xlim=c(75000, 750000)) + # Reduce horizontal whitespace around UK
theme_void() + 
theme(legend.title=element_text(size=14, vjust=1), legend.text=element_text(size=12))
ggsave("st006_03_graphical_abstract_cad.png", p2, width=7.5, height=7.5, dpi=300, units="cm")


#----
# Create merged figure
#----

g1 <- ggarrange(plotlist=plots, ncol=2, nrow=2, labels="auto", font.label=list(size=14))
ggsave("st006_01_figure1.png", g1, width=17, height=15, dpi=1200, units="cm")
system("convert st006_01_figure1.png -format tif st006_01_figure1.tif")