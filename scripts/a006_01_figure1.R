#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------------------------------------
# Setup environment
#----

options(width=200)
options(digits=6)
library(data.table)
library(ggplot2)
library(ggpubr)
library(viridis)

#----
# Load data
#----

map <- fread("../p004_maps/st004_01_coordinate_data.tsv.gz")
plot_df <- fread("../p004_maps/st004_02_map_plot_data.tsv.gz")


#----
# Create individual plots
#----

plots <- list()
for (haplo in c("I1-M253", "P-M45", "E1b1b-V13", "R1b-S749")) {
	writeLines(haplo)

	# Merge data of interest with spacial coordinates on 'id' column
	sp_map_dat <- plot_df[haplogroup_short==haplo,][map,,on="id"]

	p1 <- ggplot(sp_map_dat) +
	geom_polygon(aes(x = long, y = lat, group = group, fill = pct_haplo), color=NA, size=0.1) +
	scale_fill_viridis(name=paste0(haplo,"\nprevalence")) +
	coord_fixed(xlim=c(75000, 750000)) + # Reduce horizontal whitespace around UK
	theme_void() + 
	theme(legend.title=element_text(size=14, vjust=1), legend.text=element_text(size=12))

	plots[[haplo]] <- p1
}

#----
# Create merged figure
#----

g1 <- ggarrange(plotlist=plots, ncol=2, nrow=2, labels="auto", font.label=list(size=14))
ggsave("st006_01_figure1.png", g1, width=17, height=15, dpi=1200, units="cm", bg="white")
system("convert st006_01_figure1.png -format tif st006_01_figure1.tif")