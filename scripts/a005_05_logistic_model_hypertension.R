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
# Define haplogroups
#----

# Define groups to test based on minimum cases
cvd_groups <- c(ht_df[Cases >= min_cases, wilson_haplo_short], 
  htg_df[`CVD Cases` >= min_cases, wilson_haplo_short])

# Define groups to test based on minimum sample
quant_groups <- c(ht_df[Count >= min_sample, wilson_haplo_short], 
  htg_df[Count >= min_sample, wilson_haplo_short])


#----
# Hypertension
#----

# Convert to TRUE/FALSE
for (trait in c("hypertension","fath_hypertension", "moth_hypertension", "bro_hypertension")) {
  ph[,trait] <- ifelse(ph[,trait] > 1, TRUE, FALSE)
}

# Subset data frame
y <- c("hypertension", "fath_hypertension", "moth_hypertension", "bro_hypertension") # Dependent variables
covars <- c("haplogroup_short","array","age","fath_age","moth_age",geovars,paste0("pc",1:n_pcs)) # Haplogroup and covariates
ph1 <- ph[complete.cases(ph[,covars]),c("iid",y,covars)] # Subset to complete covariates
fwrite(ph1, "st005_05_hypertension_geo_vars.tsv", sep="\t", quote=TRUE, na="NA") # Write to file for later retrieval


# Calculate P threshold
#----------------------

# Create haplogroup variabels
for(haplo in quant_groups) {
  if(haplo %in% htgrouping$Group) {
    ph1[,haplo] <- ph1$haplogroup_short %in% htgrouping[Group==haplo, Members]
  } else {
    ph1[,haplo] <- ph1$haplogroup_short == haplo
  }
}

# Get number of independent haplogroups
ph2 <- ph1[,quant_groups]
pca <- princomp(cor(ph2, use="pair"))
cumulative_variance_explained <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
n_independent_haplos <- min(which(cumulative_variance_explained > 0.95))


#----
# Base model
#----

# Run model
#----------

res_no_geovars <- data.table()

for(haplo in quant_groups) {
  writeLines(haplo)
  fs <- paste0("hypertension ~ `",haplo,"` + array + age + I(age^2) + ",paste0("pc",1:n_pcs, collapse=" + "))
  ms <- speedglm(as.formula(fs), family=binomial(link="logit"), data=ph1, model=TRUE)
  ss <- summary(ms)
  son_coef <- ss$coefficients[grepl(paste0("^`",haplo),rownames(ss$coefficients)),,drop=FALSE]
  names(son_coef) <- c("beta","se","z","p")
  
  res1 <- data.table(
      coef=gsub("`","",gsub("TRUE","",rownames(son_coef))),
      n_haplo=sum(ms$model[,haplo]),
      n_cases=sum(ms$model[,"hypertension"]),
      n_controls=sum(!ms$model[,"hypertension"]),
      n_total=nrow(ms$model), son_coef,
      par="son"
  )

  res_no_geovars <-  rbind(res_no_geovars, res1)
}

# Validate significant
#---------------------

sig_haplos <- res_no_geovars[p < 0.1, coef]

for(haplo in sig_haplos) {
  writeLines(haplo)
  ff <- paste0("fath_hypertension ~ `",haplo,"` + array + fath_age + I(fath_age^2) + ",paste0("pc",1:n_pcs, collapse=" + "))
  mf <- speedglm(as.formula(ff), family=binomial(link="logit"), data=ph1, model=TRUE)
  sf <- summary(mf)
  fath_coef <- sf$coefficients[grepl(paste0("^`",haplo),rownames(sf$coefficients)),,drop=FALSE]
  names(fath_coef) <- c("beta","se","z","p")

  fm <- paste0("moth_hypertension ~ `",haplo,"` + array + moth_age + I(moth_age^2) + ",paste0("pc",1:n_pcs, collapse=" + "))
  mm <- speedglm(as.formula(fm), family=binomial(link="logit"), data=ph1, model=TRUE)
  sm <- summary(mm)
  moth_coef <- sm$coefficients[grepl(paste0("^`",haplo),rownames(sm$coefficients)),,drop=FALSE]
  names(moth_coef) <- c("beta","se","z","p")

  res1f <- data.table(
      coef=gsub("`","",gsub("TRUE","",rownames(fath_coef))),
      n_haplo=sum(mf$model[,haplo]),
      n_cases=sum(mf$model[,"fath_hypertension"]),
      n_controls=sum(!mf$model[,"fath_hypertension"]),
      n_total=nrow(mf$model), fath_coef,
      par="fath"
  )

  res1m <- data.table(
      coef=gsub("`","",gsub("TRUE","",rownames(moth_coef))),
      n_haplo=sum(mm$model[,haplo]),
      n_cases=sum(mm$model[,"moth_hypertension"]),
      n_controls=sum(!mm$model[,"moth_hypertension"]),
      n_total=nrow(mm$model), moth_coef,
      par="moth"
  )

  res_no_geovars <-  rbind(res_no_geovars,res1f,res1m)
}


# Adjust for multiple testing
res_no_geovars[,p_adj:=pmin(1, p * n_independent_haplos)]
fwrite(res_no_geovars, "st005_05_hypertension_base_model.csv", sep=",", na="NA", quote=FALSE)


#----
# Geovars model
#----

# Run model
#----------

res <- data.table()

for(haplo in quant_groups) {
  writeLines(haplo)
  fs <- paste0("hypertension ~ `",haplo,"` + array + age + I(age^2) + ",paste0(geovars,collapse=" + ")," + ",paste0("pc",1:n_pcs, collapse=" + "))
  ms <- speedglm(as.formula(fs), family=binomial(link="logit"), data=ph1, model=TRUE)
  ss <- summary(ms)
  son_coef <- ss$coefficients[grepl(paste0("^`",haplo),rownames(ss$coefficients)),,drop=FALSE]
  names(son_coef) <- c("beta","se","z","p")
  
  res1 <- data.table(
      coef=gsub("`","",gsub("TRUE","",rownames(son_coef))),
      n_haplo=sum(ms$model[,haplo]),
      n_cases=sum(ms$model[,"hypertension"]),
      n_controls=sum(!ms$model[,"hypertension"]),
      n_total=nrow(ms$model), son_coef,
      par="son"
  )

  res <-  rbind(res, res1)
}

# Validate significant
#---------------------

sig_haplos <- res[p < 0.1, coef]

for(haplo in sig_haplos) {
  writeLines(haplo)
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

  res1f <- data.table(
      coef=gsub("`","",gsub("TRUE","",rownames(fath_coef))),
      n_haplo=sum(mf$model[,haplo]),
      n_cases=sum(mf$model[,"fath_hypertension"]),
      n_controls=sum(!mf$model[,"fath_hypertension"]),
      n_total=nrow(mf$model), fath_coef,
      par="fath"
  )

  res1m <- data.table(
      coef=gsub("`","",gsub("TRUE","",rownames(moth_coef))),
      n_haplo=sum(mm$model[,haplo]),
      n_cases=sum(mm$model[,"moth_hypertension"]),
      n_controls=sum(!mm$model[,"moth_hypertension"]),
      n_total=nrow(mm$model), moth_coef,
      par="moth"
  )

  res <-  rbind(res,res1f,res1m)
}


# Adjust for multiple testing
res[,p_adj:=pmin(1, p * n_independent_haplos)]
fwrite(res, "st005_05_hypertension_geo_model.csv", sep=",", na="NA", quote=FALSE)


#----
# Plotting
#–---

# Create geo plot
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

ggsave("st005_05_hypertension_geo_model.pdf", p1, width=6, height=3)


# Create base plot
res1 <- res_no_geovars[coef %in% res_no_geovars[p < 0.1, coef],]

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

ggsave("st005_05_hypertension_base_model.pdf", p1, width=6, height=3)


