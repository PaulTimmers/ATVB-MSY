#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------------------------------------
# Setup environment
#----

options(width=200)
options(digits=6)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(scales)

kin_names <- c(son="Subject",fath="Father",bro="Brother",moth="Mother")
cbb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#----
# Load data
#----

res <- fread("../p005_analysis/st005_04_cvd_geo_model.csv")

#----
# Create plot
#----

res1 <- res[coef %in% res[p < 0.1, coef],]

coef_levels <- res1[par=="son",][order(beta),coef]
res1[,coef:=factor(coef, levels=coef_levels)]

res1[,kin:=factor(kin_names[par], levels=kin_names)]
breaks <- res1[,pretty(c(beta + se, beta - se), n=3)]
rf <- res1[,.(coef, beta=range(breaks), kin)]

p1 <- ggplot(res1, aes(x=beta, y=coef, colour=kin)) +
geom_vline(xintercept=0, linetype=2) +
geom_errorbarh(aes(xmin=beta - qnorm(0.975)*se, xmax=beta + qnorm(0.975)*se),position=position_dodge(-0.2), height=0) + 
geom_point(aes(shape=kin),fill="white",position=position_dodge(-0.2)) + 
geom_rangeframe(data=rf, aes(x=beta, y=coef), inherit.aes=FALSE) +
scale_shape_manual(values=c(19,19,21), guide="none") +
scale_colour_manual(values=cbb_palette, name="", guide=guide_legend(override.aes = list(shape=c(19,19,21), fill=c("#000000","#E69F00","#FFFFFF")) )) +
scale_x_continuous(breaks=breaks) +
labs(x="log(OR)", y="Haplogroup") +
theme_pubclean() + 
theme(legend.position="right")

ggsave("st006_02_figure2.png", p1, width=17, height=8.5, units="cm", dpi=1200)
system("convert st006_02_figure2.png -format tif st006_02_figure2.tif")