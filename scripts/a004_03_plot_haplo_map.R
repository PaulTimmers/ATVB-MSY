#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------------------------------------
# Setup environment
#----

options(width=200)
options(digits=6)
library(data.table)
library(ggplot2)
library(viridis)

# Variable
i <- as.numeric(commandArgs(T)[1])

#----
# Load data
#----

map <- fread("st004_01_coordinate_data.tsv.gz")
plot_df <- fread("st004_02_map_plot_data.tsv.gz")
common_haplos <- fread("st004_02_map_plot_haplos.txt", header=FALSE)[[1]]


#----
# Create plot
#----

haplo <- common_haplos[i]
writeLines(haplo)

# Merge data of interest with spacial coordinates on 'id' column
sp_map_dat <- plot_df[haplogroup_short==haplo,][map,,on="id"]

p1 <- ggplot(sp_map_dat) +
geom_polygon(aes(x = long, y = lat, group = group, fill = pct_haplo), color=NA, size=0.1) +
scale_fill_viridis(name="") +
coord_fixed() +
theme_void() + 
labs(title = paste(haplo, "prevalence"), subtitle = "By area of birth")

# Save
ggsave(paste0("st004_03_",gsub(".","_",make.names(sub("*","S",haplo,fixed=TRUE)),fixed=TRUE),"_birth.png"), p1, width=7, height=9, dpi=300)