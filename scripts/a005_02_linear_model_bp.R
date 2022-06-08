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
kin_names <- c(son="Subject",fath="Father",moth="Mother")
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
# Blood pressure
#----

# Adjust for blood pressure medication
i_medication <- with(ph, antihypertensive_meds == 1 & !is.na(antihypertensive_meds))
ph[i_medication, "sys_blood_pressure"] <- ph[i_medication, "sys_blood_pressure"] + 15 
ph[i_medication, "dia_blood_pressure"] <- ph[i_medication, "dia_blood_pressure"] + 10


# Define midpoint phenotype
ph$med_blood_pressure <- ph$dia_blood_pressure + (ph$sys_blood_pressure+ph$dia_blood_pressure)/3


# Subset data frame
y <- paste0(c("sys","dia","med"),"_blood_pressure") # Dependent variables
covars <- c("haplogroup_short","array","age",geovars,paste0("pc",1:n_pcs)) # Haplogroup and covariates
ph1 <- ph[complete.cases(ph[,covars]),c("iid",y,covars)] # Subset to complete covariates
fwrite(ph1, "st005_02_bp_geo_vars.tsv", sep="\t", quote=TRUE, na="NA") # Write to file for later retrieval


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

# Get number of independent y variables
ph3 <- ph1[,y]
pca <- princomp(cor(ph3, use="pair"))
cumulative_variance_explained <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
n_independent_vars <- min(which(cumulative_variance_explained > 0.95))

# Calculate number of independent tests
n_independent_tests <- n_independent_vars * n_independent_haplos

# Write to file
t1 <- data.table(trait="bp", n_independent_haplos,n_independent_vars,n_independent_tests)
fwrite(t1, "st005_02_bp_geo_ntests.csv", sep=",", na="NA", quote=FALSE)


#----
# Base model
#----

# Run model
#----------

res_no_geovars <- data.table()
for(haplo in quant_groups) {
  writeLines(haplo)
  fs <- paste0("sys_blood_pressure ~ `",haplo,"` + array + age + I(age^2) + ",paste0("pc",1:n_pcs, collapse=" + "))
  ms <- speedlm(as.formula(fs), data=ph1, model=TRUE)
  ss <- summary(ms)
  sys_coef <- ss$coefficients[grepl(paste0("^`",haplo),rownames(ss$coefficients)),,drop=FALSE]
  
  fd <- paste0("dia_blood_pressure ~ `",haplo,"` + array + age + I(age^2) + ",paste0("pc",1:n_pcs, collapse=" + "))
  md <- speedlm(as.formula(fd), data=ph1, model=TRUE)
  sd <- summary(md)
  dia_coef <- sd$coefficients[grepl(paste0("^`",haplo),rownames(sd$coefficients)),,drop=FALSE]

  fm <- paste0("med_blood_pressure ~ `",haplo,"` + array + age + I(age^2) + ",paste0("pc",1:n_pcs, collapse=" + "))
  mm <- speedlm(as.formula(fm), data=ph1, model=TRUE)
  sm <- summary(mm)
  med_coef <- sm$coefficients[grepl(paste0("^`",haplo),rownames(sm$coefficients)),,drop=FALSE]

  res_no_geovars <-  rbind(res_no_geovars,
  data.table(coef=rownames(sys_coef), n_haplo=sum(ms$model[,haplo]), n_total=ss$nobs, sys_coef, trait="SBP"),
  data.table(coef=rownames(dia_coef), n_haplo=sum(md$model[,haplo]), n_total=sd$nobs, dia_coef, trait="DBP"),
  data.table(coef=rownames(med_coef), n_haplo=sum(mm$model[,haplo]), n_total=sm$nobs, med_coef, trait="MAP")
  )
}

names(res_no_geovars)[c(4,7)] <- c("beta","p")
res_no_geovars[,coef:=gsub("`","",gsub("TRUE","",coef))]


# Adjust for multiple testing
res_no_geovars[,p_adj:=pmin(1, p * n_independent_tests)]
fwrite(res_no_geovars, "st005_02_bp_base_model.csv", sep=",", na="NA", quote=FALSE)




#----
# Geovars model
#----

# Run model
#----------

res <- data.table()
for(haplo in quant_groups) {
  writeLines(haplo)
  fs <- paste0("sys_blood_pressure ~ `",haplo,"` + array + age + I(age^2) + ",paste0(geovars,collapse=" + ")," + ",paste0("pc",1:n_pcs, collapse=" + "))
  ms <- speedlm(as.formula(fs), data=ph1, model=TRUE)
  ss <- summary(ms)
  sys_coef <- ss$coefficients[grepl(paste0("^`",haplo),rownames(ss$coefficients)),,drop=FALSE]
  
  fd <- paste0("dia_blood_pressure ~ `",haplo,"` + array + age + I(age^2) + ",paste0(geovars,collapse=" + ")," + ",paste0("pc",1:n_pcs, collapse=" + "))
  md <- speedlm(as.formula(fd), data=ph1, model=TRUE)
  sd <- summary(md)
  dia_coef <- sd$coefficients[grepl(paste0("^`",haplo),rownames(sd$coefficients)),,drop=FALSE]

  fm <- paste0("med_blood_pressure ~ `",haplo,"` + array + age + I(age^2) + ",paste0(geovars,collapse=" + ")," + ",paste0("pc",1:n_pcs, collapse=" + "))
  mm <- speedlm(as.formula(fm), data=ph1, model=TRUE)
  sm <- summary(mm)
  med_coef <- sm$coefficients[grepl(paste0("^`",haplo),rownames(sm$coefficients)),,drop=FALSE]

  res <-  rbind(res,
  data.table(coef=rownames(sys_coef), n_haplo=sum(ms$model[,haplo]), n_total=ss$nobs, sys_coef, trait="SBP"),
  data.table(coef=rownames(dia_coef), n_haplo=sum(md$model[,haplo]), n_total=sd$nobs, dia_coef, trait="DBP"),
  data.table(coef=rownames(med_coef), n_haplo=sum(mm$model[,haplo]), n_total=sm$nobs, med_coef, trait="MAP")
  )
}

names(res)[c(4,7)] <- c("beta","p")
res[,coef:=gsub("`","",gsub("TRUE","",coef))]


# Adjust for multiple testing
res[,p_adj:=pmin(1, p * n_independent_tests)]
fwrite(res, "st005_02_bp_geo_model.csv", sep=",", na="NA", quote=FALSE)



#----
# Plotting
#----

res1 <- res[coef %in% res[p < 0.1, coef]]
res1[,score:=sum(beta),by="coef"]
res1 <- res1[order(score)]
res1[,coef:=factor(coef, levels=unique(coef))]
breaks <- res1[,pretty(c(beta + se, beta - se), n=3)]
rf <- res1[,.(coef, beta=range(breaks), trait)]


p1 <- ggplot(res1, aes(x=beta, y=coef, colour=trait)) +
geom_vline(xintercept=0, linetype=2) +
geom_errorbarh(aes(xmin=beta - qnorm(0.975)*se, xmax=beta + qnorm(0.975)*se),position=position_dodge(-0.25), height=0) + 
geom_point(position=position_dodge(-0.25)) + 
geom_text(data=res1[p_adj < 0.05,.(beta=mean(beta),trait),by="coef"], label="*", colour="black", nudge_x=1, size=4) +
geom_rangeframe(data=rf, aes(x=beta, y=coef), inherit.aes=FALSE) +
scale_colour_manual(values=rev(cbb_palette), name="") +
scale_x_continuous(breaks=breaks) +
labs(x="mmHg", y="Haplogroup") +
theme_pubclean() + 
theme(legend.position="right")

ggsave("st005_02_bp_geo_model.pdf", p1, width=6, height=4)


