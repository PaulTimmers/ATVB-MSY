#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------------------------------------
# Setup environment
#----

source("../scripts/print_script_name.R")
options(width=200)
options(digits=6)

# Data management
library(data.table)
library(readxl)

# Stats
library(survival)
library(speedglm)
library(meta)

# Plotting
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(scales)


# Static variables
n_pcs <- 40
min_cases <- 40
min_sample <- 100

cbb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
kin_names <- c(son="Subject",fath="Father",bro="Brother",moth="Mother")
geovars <- c("northing_assessment","easting_assessment","townsend_lsoa_assessment","northing_birth","easting_birth","townsend_lsoa_birth")


# Functions
heading <- function(sentence) {
  writeLines(paste("\n\n\n=======================\n\n",sentence,"\n==========================="))
}

pvalue.extreme <- function(z) {
   log.pvalue <- log(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)
   log10.pvalue <- log.pvalue/log(10)
   mantissa <- 10^(log10.pvalue %% 1)
   exponent <- log10.pvalue %/% 1
   return(sprintf("%.2fe%03d",mantissa,exponent))
}


#----------------------------------------------------------------------------------------------------------------
# START
#----

sessionInfo()

#----
# Load data
#----

# Load UK Biobank data
base <- fread("../p003_phenotypes/st003_01_ukbb_base_data.tsv.gz")
geography <- fread("../p003_phenotypes/st003_02_ukbb_geo_data.tsv.gz")
biochem <- fread("../p003_phenotypes/st003_03_ukbb_biochem_data.tsv.gz")
haplos <- fread("../p003_phenotypes/st003_04_ukbb_haplogroups.tsv.gz")

# Merge
ph <- geography[base,,on="iid"]
ph <- biochem[ph,,on="iid"]
ph <- haplos[ph,,on="iid"]
ph <- data.frame(ph, check.names=FALSE)


# Load haplogroup counts
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
# IJK
#----

haplo <- "IJK-S137"

for (trait in c("hypertension","fath_hypertension", "moth_hypertension", "bro_hypertension")) {
  ph[,trait] <- ifelse(ph[,trait] > 1, TRUE, FALSE)
}

# Subset data frame
y <- c("hypertension", "fath_hypertension", "moth_hypertension", "bro_hypertension") # Dependent variables
covars <- c("haplogroup_short","array","age","fath_age","moth_age",geovars,paste0("pc",1:n_pcs)) # Haplogroup and covariates
ph1 <- ph[complete.cases(ph[,covars]),c("iid",y,covars)] # Subset to complete covariates


if(haplo %in% htgrouping$Group) {
  ph1[,haplo] <- ph1$haplogroup_short %in% htgrouping[Group==haplo, Members]
} else {
  ph1[,haplo] <- ph1$haplogroup_short == haplo
}

fs <- paste0("hypertension ~ `",haplo,"` + array + age + I(age^2) + ",paste0(geovars,collapse=" + ")," + ",paste0("pc",1:n_pcs, collapse=" + "))
ms <- speedglm(as.formula(fs), family=binomial(link="logit"), data=ph1, model=TRUE)
ss <- summary(ms)
son_coef <- ss$coefficients[grepl(paste0("^`",haplo),rownames(ss$coefficients)),,drop=FALSE]
names(son_coef) <- c("beta","se","z","p")

ff <- paste0("fath_hypertension ~ `",haplo,"` + array + fath_age + I(fath_age^2) + ",paste0(geovars,collapse=" + ")," + ",paste0("pc",1:n_pcs, collapse=" + "))
mf <- speedglm(as.formula(ff), family=binomial(link="logit"), data=ph1, model=TRUE)
sf <- summary(mf)
fath_coef <- sf$coefficients[grepl(paste0("^`",haplo),rownames(sf$coefficients)),,drop=FALSE]
names(fath_coef) <- c("beta","se","z","p")

fm <- paste0("moth_hypertension ~ `",haplo,"` + array + moth_age + I(moth_age^2) + ",paste0(geovars,collapse=" + ")," + ",paste0("pc",1:n_pcs, collapse=" + "))
mm <- speedglm(as.formula(fm), family=binomial(link="logit"), data=ph1, model=TRUE)
sm <- summary(mm)
moth_coef <- sm$coefficients[grepl(paste0("^`",haplo),rownames(sm$coefficients)),,drop=FALSE]
names(moth_coef) <- c("beta","se","z","p")

fb <- paste0("bro_hypertension ~ `",haplo,"` + array + age + I(age^2) + ",paste0(geovars,collapse=" + ")," + ",paste0("pc",1:n_pcs, collapse=" + "))
mb <- speedglm(as.formula(fb), family=binomial(link="logit"), data=ph1, model=TRUE)
sb <- summary(mb)
bro_coef <- sb$coefficients[grepl(paste0("^`",haplo),rownames(sb$coefficients)),,drop=FALSE]
names(bro_coef) <- c("beta","se","z","p")

res <-  rbind(
  data.table(coef=rownames(son_coef), n_haplo=sum(ms$model[,haplo]), n_cases=sum(ms$model[,"hypertension"]), n_controls=sum(!ms$model[,"hypertension"]), n_total=nrow(ms$model), son_coef, par="son"),
  data.table(coef=rownames(fath_coef), n_haplo=sum(mf$model[,haplo]), n_cases=sum(mf$model[,"fath_hypertension"]), n_controls=sum(!mf$model[,"fath_hypertension"]), n_total=nrow(mf$model), fath_coef, par="fath"),
  data.table(coef=rownames(moth_coef), n_haplo=sum(mm$model[,haplo]), n_cases=sum(mm$model[,"moth_hypertension"]), n_controls=sum(!mm$model[,"moth_hypertension"]), n_total=nrow(mm$model), moth_coef, par="moth"),
  data.table(coef=rownames(bro_coef), n_haplo=sum(mb$model[,haplo]), n_cases=sum(mb$model[,"bro_hypertension"]), n_controls=sum(!mb$model[,"bro_hypertension"]), n_total=nrow(mb$model), bro_coef, par="bro")
)
res[,coef:=gsub("`","",gsub("TRUE","",coef))]

fwrite(res, "st005_05_hypertension_geo_model_ijk.csv", sep=",", quote=FALSE, na="NA")

# IJK plot (no brother)
#----------------------

res1 <- res[par %in% c("son","fath","moth"),]
res1[,kin:=factor(kin_names[par], levels=kin_names)]
breaks <- res1[,pretty(c(beta + qnorm(0.975)*se, beta - qnorm(0.975)*se), n=3)]
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

ggsave("st005_05_hypertension_geo_model_ijk1.pdf", p1, width=6, height=2)


# IJK plot (brother)
#-------------------

res1 <- copy(res)
res1[,kin:=factor(kin_names[par], levels=kin_names)]
breaks <- res1[,pretty(c(beta + qnorm(0.975)*se, beta - qnorm(0.975)*se), n=3)]
rf <- res1[,.(coef, beta=range(breaks), kin)]

p1 <- ggplot(res1, aes(x=beta, y=coef, colour=kin)) +
geom_vline(xintercept=0, linetype=2) +
geom_errorbarh(aes(xmin=beta - qnorm(0.975)*se, xmax=beta + qnorm(0.975)*se),position=position_dodge(-0.2), height=0) + 
geom_point(aes(shape=kin),fill="white",position=position_dodge(-0.2)) + 
geom_rangeframe(data=rf, aes(x=beta, y=coef), inherit.aes=FALSE) +
scale_shape_manual(values=c(19,19,19,21), guide="none") +
scale_colour_manual(values=cbb_palette[c(1,2,4,3)], name="", guide=guide_legend(override.aes = list(shape=c(19,19,19,21), fill=c(cbb_palette[c(1,2,4)],"#FFFFFF")) )) +
scale_x_continuous(breaks=breaks) +
labs(x="log(OR)", y="Haplogroup") +
theme_pubclean() + 
theme(legend.position="right")

ggsave("st005_05_hypertension_geo_model_ijk2.pdf", p1, width=6, height=2)
 