#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------------------------------------
# Setup environment
#----

options(width=200)
options(digits=6)

# Map manipulation
library(data.table)
library(rgdal)
library(sp)
library(rgeos)
library(broom)

sessionInfo()

#----------------------------------------------------------------------------------------------------------------
# START
#----

# Load regions
#-------------

# Load data
writeLines("\n\nLoading raw data...")
sp_ward <- readOGR("../data/d003_ukbb_geography/","infuse_ward_lyr_2011_clipped") # Load Electoral Wards

# Annotate Electoral Wards with UK Country
region_mapper <- fread("../data/d003_ukbb_geography/all-areas-lookup.csv")[,.(geo_code=WED_GSS_ID, country=CTRY_NAME)]
region_mapper <- region_mapper[!duplicated(paste0(geo_code, country)), ]

# Subset to Great Britain only
sp_ward@data$country <- region_mapper[match(sp_ward@data$geo_code, geo_code), country]
sp_ward <- sp_ward[sp_ward$country %in% c("England", "Wales", "Scotland"),]

# Simplify complex coastline edges
writeLines("Simplifying coastline edges...")
sp_ward <- SpatialPolygonsDataFrame(gSimplify(sp_ward, tol=0.001, topologyPreserve=TRUE), data=sp_ward@data)

# Ignore self-intersection topology errors (if they exist)
writeLines("Resolving self-intersection topology errors...")
sp_ward <- rgeos::gBuffer(sp_ward, byid=TRUE, width=0) 

# Reformat for plotting
writeLines("\n\nReformatting spatial data...")
sp_ward_map <- data.table(broom::tidy(sp_ward, region="geo_code"))

# Write to file
fwrite(sp_ward_map, "st004_01_coordinate_data.tsv.gz", compress="gzip", sep="\t", na="NA", quote=FALSE)

writeLines("\n\n--- DONE ---\n")