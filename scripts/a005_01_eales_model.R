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
cbb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


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

count <- function(x, na.rm=T) sum(!is.na(x))
missing <- function(x, na.rm=T) sum(is.na(x))

quant_descriptives <- function(df, quant_vars, strat_vars) {
  desc_all <- data.table()
  dt <- data.table(df)[,.SD,.SDcols=c(strat_vars,quant_vars)]
  
  # Coerce classes
  set(dt, j=strat_vars, value=dt[,lapply(.SD, as.factor), .SDcols=strat_vars])
  set(dt, j=quant_vars, value=dt[,lapply(.SD, as.numeric), .SDcols=quant_vars])

  # Calculate descriptives
  for(fun in c("count","missing","mean","median","IQR")) {
    FUN <- get(fun)  
    desc <- dt[,lapply(.SD, FUN, na.rm=T),.SDcols=quant_vars]
    desc <- cbind(matrix("ALL", nrow=1, ncol=length(strat_vars), dimnames=list(fun, strat_vars)), desc)
    for(i in 1:length(strat_vars)) {
      desc1 <- dt[,lapply(.SD, FUN, na.rm=T),by=eval(strat_vars[i]),.SDcols=quant_vars]
      desc1 <- cbind(matrix("ALL", nrow=1, ncol=length(strat_vars[-i]), dimnames=list(fun, strat_vars[-i])), desc1)
      desc <- rbind(desc, desc1, fill=T)
    }
    desc[,stat:=factor(fun, levels=fun)]
    desc_all <- rbind(desc_all, desc)
  }
  return(desc_all)
}



factor_descriptives <- function(df, factor_vars, strat_vars) {
  desc_all <- data.table()
  dt <- data.table(df)[,.SD,.SDcols=c(strat_vars,factor_vars)]
  
  # Coerce classes
  set(dt, j=c(strat_vars,factor_vars), value=dt[,lapply(.SD, as.factor), .SDcols=c(strat_vars,factor_vars)])

  # Calculate descriptives
  for(factor_var in factor_vars) {
    desc <- data.frame(trait=factor_var, table(ph_f1[,c(factor_var)], useNA="always"))
    desc <- cbind(matrix("ALL", nrow=1, ncol=length(strat_vars), dimnames=list(1, strat_vars)), desc)
    names(desc)[3:4] <- c("value","count")

    desc1 <- data.frame(trait=factor_var, t(table(ph_f1[,c(factor_var,strat_vars)], useNA="always")))
    names(desc1)[3:4] <- c("value","count")

    desc <- rbind(desc, desc1)
    desc_all <- rbind(desc, desc_all)
  }
  
  return(desc_all)
}



#----------------------------------------------------------------------------------------------------------------
# START
#----

sessionInfo()

#----
# Load data
#----

base <- fread("../p003_phenotypes/st003_01_ukbb_base_data.tsv.gz")
haplos <- fread("../p003_phenotypes/st003_04_ukbb_haplogroups.tsv.gz")
ph <- data.frame(haplos[base,,on="iid"], check.names=FALSE)



#----
# Format Eales
#----

eales <- fread("../data/d004_i1_model_coefficients.csv")

# Convert to beta scale
eales[,beta:=log(or)] 
eales[,lci:=log(as.numeric(sub("-.*$","",ci)))] 
eales[,uci:=log(as.numeric(sub("^.*-","",ci)))]

# Infer standard errors
eales[p>0,z:=qnorm(p/2)]
eales[,se:=abs(beta/z)]
eales[se==0 | is.na(se),se:=pmax(uci-beta,beta-lci)/qnorm(0.975)]

# Tidy up
eales <- eales[,.(variable, beta, se, p, study="Eales 2019")]


#----
# I1 Haplogroup
#----

htgrouping <- data.table(read_xlsx("../p003_phenotypes/st003_04_haplo_groups.xlsx", sheet="Group Members"))
i1_haplos <- htgrouping[Group=="I1-M253",]
ph$i1 <- ph$haplogroup_short %in% i1_haplos$Members

writeLines("\n\nThe number of unrelated UK Biobank men carrying the I1-M253 haplogroup:")
print(with(ph, table(i1, useNA="always")))


# Identical model
#–---------------

heading("Eales et al. identical model")

# Fit model
ph$array <- factor(ph$array, levels=c("UKBL","UKBB"))
f1 <- paste0("cvd ~ i1 + age + bmi + array + I(as.factor(hypertension)) + I(as.factor(n_days_moderate_exercise)) + I(as.factor(avg_houshold_income)) + I(as.factor(smoking_history)) + I(as.factor(further_education)) + I(as.factor(current_employment)) + I(as.factor(alcohol_intake_freq)) + I(as.factor(fath_heart_disease)) + I(as.factor(moth_heart_disease)) + ", paste0("pc",1:5,collapse=" + "))
m1 <- speedglm(formula=as.formula(f1), data=ph, family=binomial(link="logit"))
s1 <- summary(m1)
a <- s1$coefficients[-1,]

print(s1)


# Write phenotype descriptives to file
vars <- all.vars(as.formula(f1))
factor_vars <- vars[grepl("array|hypertension|exercise|income|smoking|education|employment|alcohol|disease", vars)]

ph_f1 <- ph[,vars]
desc_quant <- quant_descriptives(ph_f1, quant_vars=c("age","bmi"), strat_vars="cvd")

desc_table <- rbind(
  cbind(trait="age", dcast(desc_quant, stat~cvd, value.var="age")),
  cbind(trait="bmi",dcast(desc_quant, stat~cvd, value.var="bmi")))
fwrite(desc_table, "st005_01_eales_model_descriptives.1.csv", sep=",", quote=FALSE, na="missing")

desc_factor <- factor_descriptives(ph_f1, factor_vars=factor_vars, strat_vars="cvd")
desc_table2 <- dcast(desc_factor, trait+value~cvd, value.var="count")
fwrite(desc_table2, "st005_01_eales_model_descriptives.2.csv", sep=",", quote=FALSE, na="missing")


# Get number of cases, controls, and exclusions
x <- model.frame(as.formula(f1), data=ph)
cvd_sample <- data.frame(n_cases=sum(x$cvd, na.rm=T), n_controls=sum(x$cvd==FALSE, na.rm=T), n_exclusions=nrow(ph)-nrow(x))
writeLines("\nThe CAD sample used in the identifical Eales model:")
print(format(cvd_sample,big.mark=","), row.names=FALSE)

# Format results to match Eales et al.
a$variable <- eales[,variable]
a <- data.table(a)[,.(variable, beta=as.numeric(as.character(Estimate)), se=`Std. Error`, p=pvalue.extreme(Estimate/`Std. Error`), study="This study")]


# Write effect estimate
writeLines("\nThe effect estimate of I1-M253:")
v1 <- c(a[1, exp(beta + c(0, -1, 1) * qnorm(0.975) * se)])
s1 <- "%0.2f (95%% CI %0.2f-%0.2f"
cat(do.call(sprintf, c(s1, as.list(v1))))
cat(paste0("; P = ",a[1,p],")\n"))





# Principal component sensitivity
#----------------

heading("Eales et al. model with 40 PCs")

# Fit model with 40 PCs
f2 <- paste0("cvd ~ i1 + age + bmi + array + hypertension + I(as.factor(n_days_moderate_exercise)) + I(as.factor(avg_houshold_income)) + I(as.factor(smoking_history)) + I(as.factor(further_education)) + I(as.factor(current_employment)) + I(as.factor(alcohol_intake_freq)) + fath_heart_disease + moth_heart_disease + ", paste0("pc",1:40,collapse=" + "))
m2 <- speedglm(formula=as.formula(f2), data=ph, family=binomial(link="logit"))
s2 <- summary(m2)
a2 <- s2$coefficients[-1,]

print(s2)

# Get number of cases, controls, and exclusions
x <- model.frame(as.formula(f2), data=ph)
cvd_sample <- data.frame(n_cases=sum(x$cvd, na.rm=T), n_controls=sum(x$cvd==FALSE, na.rm=T), n_exclusions=nrow(ph)-nrow(x))
writeLines("\nThe CAD sample used in the Eales model with 40 PCs:")
print(format(cvd_sample,big.mark=","), row.names=FALSE)

# Format results to match Eales et al.
a2$variable <- c(eales[,variable],paste0("PC",6:40))
a2 <- data.table(a2)[,.(variable, beta=as.numeric(as.character(Estimate)), se=`Std. Error`, p=pvalue.extreme(Estimate/`Std. Error`), study="This study")]

# Write effect estimate
writeLines("\nThe effect estimate of I1-M253:")
v2 <- c(a2[1, exp(beta + c(0, -1, 1) * qnorm(0.975) * se)])
s2 <- "%0.2f (95%% CI %0.2f-%0.2f"
cat(do.call(sprintf, c(s2, as.list(v2))))
cat(paste0("; P = ",a2[1,p],")\n"))

writeLines("\nSignificant PCs:")
v2pc <- a2[grepl("^PC", variable) & as.numeric(p) < 0.05, .(
  variable,
  or=exp(beta),
  lci= exp(beta - qnorm(0.975) * se),
  uci= exp(beta + qnorm(0.975) * se),
  p
  )]

print(v2pc)


# Parental heart disease sensitivity
#----------------

heading("Eales et al. model no father heart disease")

# Fit model without parental heart disease
f3 <- paste0("cvd ~ i1 + age + bmi + array + hypertension + I(as.factor(n_days_moderate_exercise)) + I(as.factor(avg_houshold_income)) + I(as.factor(smoking_history)) + I(as.factor(further_education)) + I(as.factor(current_employment)) + I(as.factor(alcohol_intake_freq)) + moth_heart_disease +", paste0("pc",1:40,collapse=" + "))
m3 <- speedglm(formula=as.formula(f3), data=ph, family=binomial(link="logit"))
s3 <- summary(m3)
a3 <- s3$coefficients[-1,]

print(s3)

# Get number of cases, controls, and exclusions
x <- model.frame(as.formula(f3), data=ph)
cvd_sample <- data.frame(n_cases=sum(x$cvd, na.rm=T), n_controls=sum(x$cvd==FALSE, na.rm=T), n_exclusions=nrow(ph)-nrow(x))
writeLines("\nThe CAD sample used in the Eales model without father heart disease:")
print(format(cvd_sample,big.mark=","), row.names=FALSE)

# Format results to match Eales et al.
a3$variable <- c(eales[!grepl("Paternal",a$variable),variable],paste0("PC",6:40))
a3 <- data.table(a3)[,.(variable, beta=as.numeric(as.character(Estimate)), se=`Std. Error`, p=pvalue.extreme(Estimate/`Std. Error`), study="This study")]

# Write effect estimate
writeLines("\nThe effect estimate of I1-M253:")
v3 <- c(a3[1, exp(beta + c(0, -1, 1) * qnorm(0.975) * se)])
s3 <- "%0.2f (95%% CI %0.2f-%0.2f"
cat(do.call(sprintf, c(s3, as.list(v3))))
cat(paste0("; P = ",a3[1,p],")\n"))



#----
# Plot comparison
#----

model_dat <- rbind(a,eales)
model_dat[,variable:=factor(variable, levels=rev(unique(variable)))]
rf <- model_dat[,.(variable, beta=c(-1,1))]

p1 <- ggplot(model_dat, aes(x=beta, y=variable, colour=study)) + 
geom_vline(xintercept=0, linetype=2) +
geom_errorbarh(aes(xmin=beta-qnorm(0.975)*se, xmax=beta+qnorm(0.975)*se), position=position_dodge(-0.5), alpha=0.8, height=0) +
geom_point(position=position_dodge(-0.5)) +
geom_rangeframe(data=rf, colour="black") + 
scale_colour_manual(values=cbb_palette, breaks=c("Eales 2019","This study"), name="Study") +
labs(x="log(OR)", y="") + coord_cartesian(xlim=c(-1,1)) +
theme_pubclean() +
theme(legend.position="right")

ggsave("st005_01_eales_compare.pdf", p1, width=8, height=8)


#----
# Write results to file
#----

fwrite(a, "st005_01_eales_model_identical.csv", sep=",", quote=TRUE, na="NA")
fwrite(a2, "st005_01_eales_model_40pcs.csv", sep=",", quote=TRUE, na="NA")
fwrite(a3, "st005_01_eales_model_no_fath.csv", sep=",", quote=TRUE, na="NA")

