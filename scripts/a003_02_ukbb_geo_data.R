#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------------------------------------
# Setup environment
#----

source("../scripts/print_script_name.R")
options(width=200)
options(digits=6)

# R packages

library(data.table)
library(rgdal)
library(sp)

# Variables
regions <- c("rgn", "dist", "ward", "lsoa")
region_files <- c(rgn="infuse_rgn_2011_clipped", dist="infuse_dist_lyr_2011_clipped", ward="infuse_ward_lyr_2011_clipped", lsoa="infuse_lsoa_lyr_2011_clipped")
pretty_names <- c(rgn="Region", dist="Local Authories", ward="Wards and Electoral Divisions", lsoa="Lower Layer Super Output Areas")
suffixes <- c("assessment", "birth")

traits <- c(
"f.eid",
"f.20074.0.0",
"f.20075.0.0",
"f.129.0.0",
"f.130.0.0",
"f.1647.0.0"
)

trait_names <- c(
"iid",
"easting_assessment",
"northing_assessment",
"northing_birth",
"easting_birth",
"place_birth"
)


#----------------------------------------------------------------------------------------------------------------
# Start
#----

sessionInfo()

#----
# Load data
#----

# Load coordinate phenotypes
df_phenotypes <- fread("../data/d003_ukbb_data.tsv", select=traits, col.names=trait_names)


# Load townsend regions
townsend_mapper <- fread(paste0("../data/d003_ukbb_geography/townsend_lsoa_2011.csv"), col.names=c("id","geo_code","geo_label","town_di","town_di_quintile"))


#----
# Map coordinates onto geographical areas
#----

for (region in regions) {

# Load spatial coordinates
writeLines(paste0("\n\n========\n\n",pretty_names[region],"\n===================\n"))
sp <- readOGR("../data/d003_ukbb_geography", region_files[region])

# Create geography variables
	for (suffix in suffixes) {
		
		traits <- paste(c("easting","northing"),suffix,sep="_") 
		df_phenotypes1 <- df_phenotypes[,.SD,.SDcols=c("iid",traits)] # Select traits
		names(df_phenotypes1) <- c("iid","easting","northing")

		# Select complete coordinate data
		coord_dat <- df_phenotypes1[!is.na(easting) & !is.na(northing), ]
		coord_dat_iid <- coord_dat$iid

		# Convert raw easting and northing data to spacial coordinates
		coordinates(coord_dat) <- c("easting","northing")
		proj4string(coord_dat) <- CRS("+init=epsg:27700")

		# Reproject easting and northing data to match lsoa/msoa/ward/dist coordinate system
		coord_sp <- spTransform(coord_dat, CRS(proj4string(sp)))

		# Map easting northing points onto lsoa/msoa/ward/dist polygons
		cat(paste0("Getting ",pretty_names[region]," labels from ",paste0(traits,collapse="/"),"... "))
		area_dat <- coord_sp %over% sp
		area_dat <- townsend_mapper[area_dat,,on=c("geo_code","geo_label")]
		cat("done!\n")

		# Add results onto original data frame
		df_phenotypes[match(coord_dat_iid,iid),paste(region,suffix,sep="_")] <- area_dat$geo_code
		df_phenotypes[match(coord_dat_iid,iid),paste(region,suffix,"label",sep="_")] <- area_dat$geo_label
		if(region == "lsoa") df_phenotypes[match(coord_dat_iid,iid),paste("townsend",region,suffix,sep="_")] <- area_dat$town_di

		# Write (intermediate) results
		fwrite(df_phenotypes, paste0("st003_02_ukbb_geo_data.tsv.gz"), compress="gzip", sep="\t", quote=TRUE, na="NA")
	}
}