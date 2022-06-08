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
# Lipid levels
#----

names(ph) <- gsub("cholesterol_1","cholesterol",names(ph))
names(ph) <- gsub("log_triglycerides_1","log_triglycerides",names(ph))

# Adjust for lipid-lowering medication
i_medication <- with(ph, cholesterol_meds == 1 & !is.na(cholesterol_meds))
ph[i_medication, "ldl_cholesterol"] <- ph[i_medication, "ldl_cholesterol"] / 0.7 
ph[i_medication, "total_cholesterol"] <- ph[i_medication, "total_cholesterol"] / 0.8


# Subset data frame
y <- c(paste0(c("total","ldl","hdl"),"_cholesterol"),"log_triglycerides") # Dependent variables
covars <- c("haplogroup_short","array","age",geovars,paste0("pc",1:n_pcs)) # Haplogroup and covariates
ph1 <- ph[complete.cases(ph[,covars]),c("iid",y,covars)] # Subset to complete covariates
fwrite(ph1, "st005_03_lipid_geo_vars.tsv", sep="\t", quote=TRUE, na="NA") # Write to file for later retrieval


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
t1 <- data.table(trait="lipid", n_independent_haplos,n_independent_vars,n_independent_tests)
fwrite(t1, "st005_03_lipid_geo_ntests.csv", sep=",", na="NA", quote=FALSE)

#----
# Base model
#----

# Run model
#----------

res_no_geovars <- data.table()
for(haplo in quant_groups) {
  writeLines(haplo)
  ft <- paste0("total_cholesterol ~ `",haplo,"` + array + age + I(age^2) + ",paste0("pc",1:n_pcs, collapse=" + "))
  mt <- speedlm(as.formula(ft), data=ph1, model=TRUE)
  st <- summary(mt)
  total_coef <- st$coefficients[grepl(paste0("^`",haplo),rownames(st$coefficients)),,drop=FALSE]
  
  fl <- paste0("ldl_cholesterol ~ `",haplo,"` + array + age + I(age^2) + ",paste0("pc",1:n_pcs, collapse=" + "))
  ml <- speedlm(as.formula(fl), data=ph1, model=TRUE)
  sl <- summary(ml)
  ldl_coef <- sl$coefficients[grepl(paste0("^`",haplo),rownames(sl$coefficients)),,drop=FALSE]
  
  fh <- paste0("hdl_cholesterol ~ `",haplo,"` + array + age + I(age^2) + ",paste0("pc",1:n_pcs, collapse=" + "))
  mh <- speedlm(as.formula(fh), data=ph1, model=TRUE)
  sh <- summary(mh)
  hdl_coef <- sh$coefficients[grepl(paste0("^`",haplo),rownames(sh$coefficients)),,drop=FALSE]

  fy <- paste0("log_triglycerides ~ `",haplo,"` + array + age + I(age^2) + ",paste0("pc",1:n_pcs, collapse=" + "))
  my <- speedlm(as.formula(fy), data=ph1, model=TRUE)
  sy <- summary(my)
  triglyc_coef <- sy$coefficients[grepl(paste0("^`",haplo),rownames(sy$coefficients)),,drop=FALSE]

  res_no_geovars <-  rbind(res_no_geovars,
  data.table(coef=rownames(total_coef), n_haplo=sum(mt$model[,haplo]), n_total=st$nobs, total_coef, trait="Total cholesterol"),
  data.table(coef=rownames(ldl_coef), n_haplo=sum(ml$model[,haplo]), n_total=sl$nobs, ldl_coef, trait="LDL cholesterol"),
  data.table(coef=rownames(hdl_coef), n_haplo=sum(mh$model[,haplo]), n_total=sh$nobs, hdl_coef, trait="HDL cholesterol"),
  data.table(coef=rownames(triglyc_coef), n_haplo=sum(my$model[,haplo]), n_total=sy$nobs, triglyc_coef, trait="Triglycerides")
  )
}

names(res_no_geovars)[c(4,7)] <- c("beta","p")
res_no_geovars[,coef:=gsub("`","",gsub("TRUE","",coef))]


# Adjust for multiple testing
res_no_geovars[,p_adj:=pmin(1, p * n_independent_tests)]
fwrite(res_no_geovars, "st005_03_lipid_base_model.csv", sep=",", na="NA", quote=FALSE)


#----
# Geovars model
#----

# Run model
#----------

res <- data.table()
for(haplo in quant_groups) {
  writeLines(haplo)
  ft <- paste0("total_cholesterol ~ `",haplo,"` + array + age + I(age^2) + ",paste0(geovars,collapse=" + ")," + ",paste0("pc",1:n_pcs, collapse=" + "))
  mt <- speedlm(as.formula(ft), data=ph1, model=TRUE)
  st <- summary(mt)
  total_coef <- st$coefficients[grepl(paste0("^`",haplo),rownames(st$coefficients)),,drop=FALSE]
  
  fl <- paste0("ldl_cholesterol ~ `",haplo,"` + array + age + I(age^2) + ",paste0(geovars,collapse=" + ")," + ",paste0("pc",1:n_pcs, collapse=" + "))
  ml <- speedlm(as.formula(fl), data=ph1, model=TRUE)
  sl <- summary(ml)
  ldl_coef <- sl$coefficients[grepl(paste0("^`",haplo),rownames(sl$coefficients)),,drop=FALSE]
  
  fh <- paste0("hdl_cholesterol ~ `",haplo,"` + array + age + I(age^2) + ",paste0(geovars,collapse=" + ")," + ",paste0("pc",1:n_pcs, collapse=" + "))
  mh <- speedlm(as.formula(fh), data=ph1, model=TRUE)
  sh <- summary(mh)
  hdl_coef <- sh$coefficients[grepl(paste0("^`",haplo),rownames(sh$coefficients)),,drop=FALSE]

  fy <- paste0("log_triglycerides ~ `",haplo,"` + array + age + I(age^2) + ",paste0(geovars,collapse=" + ")," + ",paste0("pc",1:n_pcs, collapse=" + "))
  my <- speedlm(as.formula(fy), data=ph1, model=TRUE)
  sy <- summary(my)
  triglyc_coef <- sy$coefficients[grepl(paste0("^`",haplo),rownames(sy$coefficients)),,drop=FALSE]

  res <-  rbind(res,
  data.table(coef=rownames(total_coef), n_haplo=sum(mt$model[,haplo]), n_total=st$nobs, total_coef, trait="Total cholesterol"),
  data.table(coef=rownames(ldl_coef), n_haplo=sum(ml$model[,haplo]), n_total=sl$nobs, ldl_coef, trait="LDL cholesterol"),
  data.table(coef=rownames(hdl_coef), n_haplo=sum(mh$model[,haplo]), n_total=sh$nobs, hdl_coef, trait="HDL cholesterol"),
  data.table(coef=rownames(triglyc_coef), n_haplo=sum(my$model[,haplo]), n_total=sy$nobs, triglyc_coef, trait="Triglycerides")
  )
}

names(res)[c(4,7)] <- c("beta","p")
res[,coef:=gsub("`","",gsub("TRUE","",coef))]


# Adjust for multiple testing
res[,p_adj:=pmin(1, p * n_independent_tests)]
fwrite(res, "st005_03_lipid_geo_model.csv", sep=",", na="NA", quote=FALSE)


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
geom_errorbarh(aes(xmin=beta - qnorm(0.975)*se, xmax=beta + qnorm(0.975)*se),position=position_dodge(-0.5), height=0) + 
geom_point(position=position_dodge(-0.5)) + 
geom_text(data=res1[p_adj < 0.05,.(beta=mean(beta),trait),by="coef"], label="*", colour="black", nudge_x=1, size=4) +
geom_rangeframe(data=rf, aes(x=beta, y=coef), inherit.aes=FALSE) +
scale_colour_manual(values=rev(cbb_palette), name="") +
scale_x_continuous(breaks=breaks) +
labs(x="Beta", y="Haplogroup") +
theme_pubclean() + 
theme(legend.position="right")

ggsave("st005_03_lipid_geo_model.pdf", p1, width=8, height=6)


