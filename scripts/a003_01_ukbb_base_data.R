#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------------------------------------
# Setup environment
#----

source("../scripts/print_script_name.R")
options(width=200)
options(digits=6)


# R packages

library(data.table)
library(reshape2)
library(lubridate)
library(stringr)


# Variables

traits <- c(
"f.eid", # iid
"f.31.0.0", # sex
"f.34.0.0", # yob
"f.52.0.0", # mob
"f.53.0.0", # date_ass_c
"f.54.0.0", # ass_c
"f.738.0.0", # income
"f.884.0.0", # moderate_phys_activity_days_per_week
"f.1558.0.0", # alcohol_intake
"f.1767.0.0", # adopted
"f.1797.0.0", # fath_alive
"f.1807.0.0", # fath_age_d
"f.1835.0.0", # moth_alive
"f.1845.0.0", # moth_age
"f.1873.0.0", # n_brothers
"f.1883.0.0", # n_sisters
"f.2946.0.0", # fath_age
"f.3526.0.0", # moth_age_d
paste0("f.4079.0.",0:1), # dia_blood_pressure
paste0("f.4080.0.",0:1), # sys_blood_pressure
"f.6138.0.0", # education
"f.6142.0.0", # employment
paste0("f.6177.0.",0:2), # cvd_meds
paste0("f.20002.0.",0:28), # ill_code
paste0("f.20003.0.",0:47), # med_code
paste0("f.20004.0.",0:31), # op_code
paste0("f.20107.0.",0:9), # fath_ill
paste0("f.20110.0.",0:10), # moth_ill
paste0("f.20111.0.",0:10), # sib_ill
"f.20116.0.0", # smoking
"f.21000.0.0", # ethnicity
"f.21001.0.0", # bmi
"f.21003.0.0", # age
"f.22000.0.0", # array
"f.22001.0.0", # genetic_sex
paste0("f.22009.0.",1:40) # pc1 - pc40

)

trait_names <- c(
"iid",
"sex",
"yob",
"mob",
"date_ass_c",
"ass_c",
"income",
"moderate_phys_activity_days_per_week",
"alcohol_intake",
"adopted",
"fath_alive",
"fath_age_d",
"moth_alive",
"moth_age",
"n_brothers",
"n_sisters",
"fath_age",
"moth_age_d",
paste0("dia_blood_pressure_",1:2),
paste0("sys_blood_pressure_",1:2),
"education",
"employment",
paste0("cvd_meds_",1:3),
paste0("ill_code_",1:29),
paste0("med_code_",1:48),
paste0("op_code_",1:32),
paste0("fath_ill_",1:10),
paste0("moth_ill_",1:11),
paste0("sib_ill_",1:11),
"smoking",
"ethnicity",
"bmi",
"age",
"array",
"genetic_sex",
paste0("pc",1:40)
)


# Functions

heading <- function(sentence) {
	writeLines(paste("\n\n\n=======================\n\n",
	sentence,"\n==========================="))
}


tukey_qc <- function(x){
  var_name <- deparse(substitute(x)) # Get variable name
  nx <- sum(!is.na(x))

  repeat {
  quant <- quantile(x, na.rm=T, probs=c(0.25,0.75)) # Get quantiles
  iqr <- diff(quant) # Interquantile range

  # Min = lower quantile minus 3x interquartile range
  # Max = upper quantile plus 3x interquartile range
  allowed <- range(quant + c(-3,3)*iqr) 

  qcx <- ifelse(x < allowed[1] | x > allowed[2], NA, x) # Remove values outside of allowed range
  if(sum(is.na(x)) == sum(is.na(qcx))) { break } else { x <- qcx } # Repeat to see if excluded values alter interquartile range
  }

  writeLines(sprintf("%s: Removed %i values outside of [%.4g to %.4g]",var_name[1],nx-sum(!is.na(qcx)),allowed[1],allowed[2]))
  return(qcx)
}




#----------------------------------------------------------------------------------------------------------------
# START
#----

sessionInfo()


#----
# Load data
#----

ph <- fread("../data/d003_ukbb_data.tsv", select=traits, col.names=trait_names, data.table=FALSE)
iids <- fread("../p001_qc_variants/st001_03_iids.txt", header=FALSE, col.names=c("iid"))


# Subset to unrelated, genomically british sample
ph <- ph[ph$iid %in% iids$iid, ]
ph <- ph[ph$sex == "Male" & ph$genetic_sex == "Male", ]



# Find adoptees
#--------------

heading("Adoptees")

ph$not_adopted <- c(`No`=TRUE,`Yes`=FALSE,`Do not know`=NA,`Prefer not to answer`=NA)[ph$adopted]

# Exclude kin trait values for adopted individuals
print(data.table(ph)[,.N,by=c("not_adopted","adopted")][order(not_adopted,-N)])

kinvars <- names(ph)[grepl("sib|fath|moth",names(ph))]
ph[!ph$not_adopted | is.na(ph$not_adopted), kinvars] <- NA



# Find full brothers
#-------------------

heading("Brothers")

writeLines(paste0("\nThere are ",nrow(ph[ph$n_brothers == 1 & ph$n_sisters ==0, ]), " individuals with a single, male sibling."))

sibvars <- names(ph)[grepl("sib",names(ph))]
ph[ph$n_brothers != 1 | is.na(ph$n_brothers) | ph$n_sisters != 0 | is.na(ph$n_sisters),sibvars] <- NA


#----------------------------------------------------------------------------------------------------------------
# Define Eales et al. covariates
#----

# Age
#----

heading("Age")
ph$age <- with(ph, tukey_qc(age))


# Array
#------

heading("Array")
ph$array <- factor(ifelse(ph$array < 0 | ph$array == 1000, "UKBL","UKBB"),levels=c("UKBL","UKBB"))
print(table(ph$array))


# BMI
#----

heading("Body mass index")
ph$bmi <- with(ph, tukey_qc(bmi))

print(summary(ph$bmi))


# Hypertension
#-------------

heading("Hypertension")

medmap <- c(`Blood pressure medication`=1,
  `Cholesterol lowering medication`=0,
  `Do not know`=NA,
  `Insulin`=0,
  `None of the above`=0,
  `Prefer not to answer`=NA
)
ph$antihypertensive_meds_1 <- medmap[ph$cvd_meds_1]
ph$antihypertensive_meds_2 <- medmap[ph$cvd_meds_2]
ph$antihypertensive_meds_3 <- medmap[ph$cvd_meds_3]

ph$antihypertensive_meds <- ph$antihypertensive_meds_1 | (ph$antihypertensive_meds_2 & !is.na(ph$antihypertensive_meds_2)) | (ph$antihypertensive_meds_3 & !is.na(ph$antihypertensive_meds_3))


writeLines("\nAntihypertensive medication:")
print(table(ph$antihypertensive_meds, useNA="ifany"))
cat("\n")
print(data.table(ph)[,.N,by=c("antihypertensive_meds_1","cvd_meds_1")][order(antihypertensive_meds_1,-N)])
print(table(ph$antihypertensive_meds_2, useNA="ifany"))
cat("\n")
print(data.table(ph)[,.N,by=c("antihypertensive_meds_2","cvd_meds_2")][order(antihypertensive_meds_2,-N)])
print(table(ph$antihypertensive_meds_3, useNA="ifany"))
cat("\n")
print(data.table(ph)[,.N,by=c("antihypertensive_meds_3","cvd_meds_3")][order(antihypertensive_meds_3,-N)])


writeLines("\nBlood pressure:")
ph$sys_blood_pressure <- with(ph, (sys_blood_pressure_1 + sys_blood_pressure_2)/2) # Take average of measurements
ph$sys_blood_pressure <- with(ph, ifelse(is.na(sys_blood_pressure), sys_blood_pressure_1, sys_blood_pressure)) # If only one measurement, take that one
ph$sys_blood_pressure <- with(ph, ifelse(is.na(sys_blood_pressure), sys_blood_pressure_2, sys_blood_pressure)) # If only one measurement, take that one
ph$dia_blood_pressure <- with(ph, (dia_blood_pressure_1 + dia_blood_pressure_2)/2) # Take average of measurements
ph$dia_blood_pressure <- with(ph, ifelse(is.na(dia_blood_pressure), dia_blood_pressure_1, dia_blood_pressure)) # If only one measurement, take that one
ph$dia_blood_pressure <- with(ph, ifelse(is.na(dia_blood_pressure), dia_blood_pressure_2, dia_blood_pressure)) # If only one measurement, take that one

ph$sys_blood_pressure <- with(ph, tukey_qc(sys_blood_pressure))
ph$dia_blood_pressure <- with(ph, tukey_qc(dia_blood_pressure))



writeLines("\nHypertension:")
ph$hypertension <- NA
ph$hypertension <- with(ph, ifelse(sys_blood_pressure >= 140 & !is.na(sys_blood_pressure), 2, hypertension))
ph$hypertension <- with(ph, ifelse(dia_blood_pressure >= 90 & !is.na(dia_blood_pressure), 2, hypertension))
ph$hypertension <- with(ph, ifelse(antihypertensive_meds == 1 & !is.na(antihypertensive_meds), 2, hypertension))
ph$hypertension <- with(ph, ifelse(sys_blood_pressure < 140 & dia_blood_pressure < 90 & antihypertensive_meds == 0, 1, hypertension))
print(table(ph$hypertension, useNA="ifany"))





# Physical activity
#------------------

heading("Physical activity")

# Number of days per week of moderate physical activity (at least 10 minutes)
ph$n_days_moderate_exercise <- with(ph, as.factor(ifelse(moderate_phys_activity_days_per_week < 0, NA, moderate_phys_activity_days_per_week+1)))

writeLines("\nNumber of days of moderate physical activity per week (coded 1-8):")
print(table(ph$n_days_moderate_exercise, useNA="ifany"))





# Household income
#-----------------

heading("Household income")

# Average total annual household income before tax.

incmap <- c(`Less than 18,000`=1,
  `18,000 to 30,999`=2,
  `31,000 to 51,999`=3,
  `52,000 to 100,000`=4,
  `Greater than 100,000`=5,
  `Do not know`=NA,
  `Prefer not to answer`=NA
  )

ph$avg_houshold_income <- as.factor(incmap[ph$income])

writeLines("\nHoushold income:")
print(table(ph$avg_houshold_income, useNA="ifany"))
cat("\n")
print(data.table(ph)[,.N,by=c("avg_houshold_income","income")][order(avg_houshold_income,-N)])



# Smoking history
#----------------

heading("Smoking")

# Positive smoking history defined as being a current smoker or being a previous smoker. Negative smoking history was defined as neither being a current or previous smoker.

smokmap <- c(Current=2, Never=1, `Prefer not to answer`=NA, Previous=2)
ph$smoking_history <- smokmap[ph$smoking]

writeLines("\nSmoking history:")
print(table(ph$smoking_history, useNA="ifany"))
cat("\n")
print(data.table(ph)[,.N,by=c("smoking_history","smoking")][order(smoking_history,-N)])



# Further education
#------------------

heading("Further education")

# Individuals who self-reported completion of A-levels/AS-level or college/university degree or other professional qualifications (e.g. nursing or teaching) were classified as having entered further education. Individuals with lower qualifications (CSEs/O-levels/GCSEs/NVQs/HNC/HND) or the absence of any formal qualifications were classified as not having entered further education.

edumap <- c(`A levels/AS levels or equivalent`=2,
            `College or University degree`=2,
            `Other professional qualifications eg: nursing, teaching`=2,
            `CSEs or equivalent`=1,
            `NVQ or HND or HNC or equivalent`=1,
            `O levels/GCSEs or equivalent`=1,
            `None of the above`=1,
            `Prefer not to answer`=NA
)

ph$further_education <- edumap[ph$education]

writeLines("\nFurther education:")
print(table(ph$further_education, useNA="ifany"))
cat("\n")
print(data.table(ph)[,.N,by=c("further_education","education")][order(further_education,-N)])





# Employment
#-----------

heading("Employment")

# Current status of employment

empmap <- c(`Doing unpaid or voluntary work`=6,
`Full or part-time student`=4,
`In paid employment or self-employed`=1,
`Looking after home and/or family`=2,
`None of the above`=NA,
`Prefer not to answer`=NA,
`Retired`=5,
`Unable to work because of sickness or disability`=3,
`Unemployed`=6
)
ph$current_employment <- as.factor(empmap[ph$employment])

writeLines("\nCurrent employement:")
print(table(ph$current_employment, useNA="ifany"))
cat("\n")
print(data.table(ph)[,.N,by=c("current_employment","employment")][order(current_employment,-N)])





# Alcohol intake frequency
#-------------------------

heading("Alcohol intake")

# Average frequency of alcohol intake

alcmap <- c(`Daily or almost daily`=6,
`Never`=1,
`Once or twice a week`=4,
`One to three times a month`=3,
`Prefer not to answer`=NA,
`Special occasions only`=2,
`Three or four times a week`=5
)

ph$alcohol_intake_freq <- as.factor(alcmap[ph$alcohol_intake])

writeLines("\nAverage frequency of alcohol intake:")
print(table(ph$alcohol_intake_freq, useNA="ifany"))
cat("\n")
print(data.table(ph)[,.N,by=c("alcohol_intake_freq","alcohol_intake")][order(alcohol_intake_freq,-N)])





# Father heart disease
#-----------------------

heading("Father heart disease")

# Father heart disease

ph$fath_heart_disease <- ph$fath_heart_disease2 <- ifelse(is.na(ph$fath_ill_1), NA, 1)
fath_ill_cols <- grep("fath_ill",names(ph),value=T)

# 1: Any reported incidence is a case; 2: Exclude NAs
for(fath_ill_col in fath_ill_cols) {
ph$fath_heart_disease <- ifelse(ph[,fath_ill_col] == "Heart disease" & !is.na(ph[,fath_ill_col]), 2, ph$fath_heart_disease)
ph$fath_heart_disease <- ifelse(ph$fath_heart_disease == 1 & ph[,fath_ill_col] %in% c("Prefer not to answer (group 1)") & !is.na(ph[,fath_ill_col]), NA, ph$fath_heart_disease)

ph$fath_heart_disease2 <- ifelse(ph[,fath_ill_col] == "Heart disease" & !is.na(ph[,fath_ill_col]), 2, ph$fath_heart_disease2)
ph$fath_heart_disease2 <- ifelse(ph$fath_heart_disease2 == 1 & ph[,fath_ill_col] %in% c("Prefer not to answer (group 1)", "Do not know (group 1)") & !is.na(ph[,fath_ill_col]), NA, ph$fath_heart_disease2)
}

writeLines("\nExcluding 'Prefer not to answer':")
print(table(ph$fath_heart_disease, useNA="ifany"))
writeLines("\nExcluding 'Prefer not to answer' and 'Do not know':")
print(table(ph$fath_heart_disease2, useNA="ifany"))


# Mother heart disease
#-----------------------

heading("Mother heart disease")

# Mother heart disease

ph$moth_heart_disease <- ph$moth_heart_disease2 <- ifelse(is.na(ph$moth_ill_1), NA, 1)
moth_ill_cols <- grep("moth_ill",names(ph),value=T)

# 1: Any reported incidence is a case; 2: Exclude NAs
for(moth_ill_col in moth_ill_cols) {
ph$moth_heart_disease <- ifelse(ph[,moth_ill_col] == "Heart disease" & !is.na(ph[,moth_ill_col]), 2, ph$moth_heart_disease)
ph$moth_heart_disease <- ifelse(ph$moth_heart_disease == 1 & ph[,moth_ill_col] %in% c("Prefer not to answer (group 1)") & !is.na(ph[,moth_ill_col]), NA, ph$moth_heart_disease)

ph$moth_heart_disease2 <- ifelse(ph[,moth_ill_col] == "Heart disease" & !is.na(ph[,moth_ill_col]), 2, ph$moth_heart_disease2)
ph$moth_heart_disease2 <- ifelse(ph$moth_heart_disease2 == 1 & ph[,moth_ill_col] %in% c("Prefer not to answer (group 1)", "Do not know (group 1)") & !is.na(ph[,moth_ill_col]), NA, ph$moth_heart_disease2)
}

writeLines("\nExcluding 'Prefer not to answer':")
print(table(ph$moth_heart_disease, useNA="ifany"))
writeLines("\nExcluding 'Prefer not to answer' and 'Do not know':")
print(table(ph$moth_heart_disease2, useNA="ifany"))


# Brother heart disease
#----------------------

heading("Brother heart disease")

# Brother heart disease

ph$bro_heart_disease <- ph$bro_heart_disease2 <- ifelse(is.na(ph$sib_ill_1), NA, 1)
sib_ill_cols <- grep("sib_ill",names(ph),value=T)

# 1: Any reported incidence is a case; 2: Exclude NAs
for(sib_ill_col in sib_ill_cols) {
ph$bro_heart_disease <- ifelse(ph[,sib_ill_col] == "Heart disease" & !is.na(ph[,sib_ill_col]), 2, ph$bro_heart_disease)
ph$bro_heart_disease <- ifelse(ph$bro_heart_disease == 1 & ph[,sib_ill_col] %in% c("Prefer not to answer (group 1)") & !is.na(ph[,sib_ill_col]), NA, ph$bro_heart_disease)

ph$bro_heart_disease2 <- ifelse(ph[,sib_ill_col] == "Heart disease" & !is.na(ph[,sib_ill_col]), 2, ph$bro_heart_disease2)
ph$bro_heart_disease2 <- ifelse(ph$bro_heart_disease2 == 1 & ph[,sib_ill_col] %in% c("Prefer not to answer (group 1)", "Do not know (group 1)") & !is.na(ph[,sib_ill_col]), NA, ph$bro_heart_disease2)
}

writeLines("\nExcluding 'Prefer not to answer':")
print(table(ph$bro_heart_disease, useNA="ifany"))
writeLines("\nExcluding 'Prefer not to answer' and 'Do not know':")
print(table(ph$bro_heart_disease2, useNA="ifany"))



# Father hypertension
#--------------------

heading("Father hypertension")

# Father hypertension

ph$fath_hypertension <- ph$fath_hypertension2 <- ifelse(is.na(ph$fath_ill_1), NA, 1)
fath_ill_cols <- grep("fath_ill",names(ph),value=T)

# 1: Any reported incidence is a case; 2: Exclude NAs
for(fath_ill_col in fath_ill_cols) {
ph$fath_hypertension <- ifelse(ph[,fath_ill_col] == "High blood pressure" & !is.na(ph[,fath_ill_col]), 2, ph$fath_hypertension)
ph$fath_hypertension <- ifelse(ph$fath_hypertension == 1 & ph[,fath_ill_col] %in% c("Prefer not to answer (group 1)") & !is.na(ph[,fath_ill_col]), NA, ph$fath_hypertension)

ph$fath_hypertension2 <- ifelse(ph[,fath_ill_col] == "High blood pressure" & !is.na(ph[,fath_ill_col]), 2, ph$fath_hypertension2)
ph$fath_hypertension2 <- ifelse(ph$fath_hypertension2 == 1 & ph[,fath_ill_col] %in% c("Prefer not to answer (group 1)", "Do not know (group 1)") & !is.na(ph[,fath_ill_col]), NA, ph$fath_hypertension2)
}

writeLines("\nExcluding 'Prefer not to answer':")
print(table(ph$fath_hypertension, useNA="ifany"))
writeLines("\nExcluding 'Prefer not to answer' and 'Do not know':")
print(table(ph$fath_hypertension2, useNA="ifany"))


# Mother hypertension
#-----------------------

heading("Mother hypertension")

# Mother hypertension

ph$moth_hypertension <- ph$moth_hypertension2 <- ifelse(is.na(ph$moth_ill_1), NA, 1)
moth_ill_cols <- grep("moth_ill",names(ph),value=T)

# 1: Any reported incidence is a case; 2: Exclude NAs
for(moth_ill_col in moth_ill_cols) {
ph$moth_hypertension <- ifelse(ph[,moth_ill_col] == "High blood pressure" & !is.na(ph[,moth_ill_col]), 2, ph$moth_hypertension)
ph$moth_hypertension <- ifelse(ph$moth_hypertension == 1 & ph[,moth_ill_col] %in% c("Prefer not to answer (group 1)") & !is.na(ph[,moth_ill_col]), NA, ph$moth_hypertension)

ph$moth_hypertension2 <- ifelse(ph[,moth_ill_col] == "High blood pressure" & !is.na(ph[,moth_ill_col]), 2, ph$moth_hypertension2)
ph$moth_hypertension2 <- ifelse(ph$moth_hypertension2 == 1 & ph[,moth_ill_col] %in% c("Prefer not to answer (group 1)", "Do not know (group 1)") & !is.na(ph[,moth_ill_col]), NA, ph$moth_hypertension2)
}

writeLines("\nExcluding 'Prefer not to answer':")
print(table(ph$moth_hypertension, useNA="ifany"))
writeLines("\nExcluding 'Prefer not to answer' and 'Do not know':")
print(table(ph$moth_hypertension2, useNA="ifany"))


# Brother hypertension
#-----------------------

heading("Brother hypertension")

# Brother hypertension

ph$bro_hypertension <- ph$bro_hypertension2 <- ifelse(is.na(ph$sib_ill_1), NA, 1)
sib_ill_cols <- grep("sib_ill",names(ph),value=T)

# 1: Any reported incidence is a case; 2: Exclude NAs
for(sib_ill_col in sib_ill_cols) {
ph$bro_hypertension <- ifelse(ph[,sib_ill_col] == "High blood pressure" & !is.na(ph[,sib_ill_col]), 2, ph$bro_hypertension)
ph$bro_hypertension <- ifelse(ph$bro_hypertension == 1 & ph[,sib_ill_col] %in% c("Prefer not to answer (group 1)") & !is.na(ph[,sib_ill_col]), NA, ph$bro_hypertension)

ph$bro_hypertension2 <- ifelse(ph[,sib_ill_col] == "High blood pressure" & !is.na(ph[,sib_ill_col]), 2, ph$bro_hypertension2)
ph$bro_hypertension2 <- ifelse(ph$bro_hypertension2 == 1 & ph[,sib_ill_col] %in% c("Prefer not to answer (group 1)", "Do not know (group 1)") & !is.na(ph[,sib_ill_col]), NA, ph$bro_hypertension2)
}

writeLines("\nExcluding 'Prefer not to answer':")
print(table(ph$bro_hypertension, useNA="ifany"))
writeLines("\nExcluding 'Prefer not to answer' and 'Do not know':")
print(table(ph$bro_hypertension2, useNA="ifany"))



#----------------------------------------------------------------------------------------------------------------
# Heart disease definition
#----

heading("Heart disease definition")


ph$cvd <- FALSE


# Self-reported
#--------------

# Myocardial infarction
ill_code_cols <- grep("ill_code",names(ph),value=T)
for(ill_code_col in ill_code_cols) {
  ph$cvd <- ifelse(ph[,ill_code_col]==1079 & !is.na(ph[,ill_code_col]), TRUE, ph$cvd)
}

# Percutaneous coronary intervention OR coronary artery bypass graft
op_code_cols <- grep("op_code",names(ph),value=T)
for(op_code_col in op_code_cols) {
  ph$cvd <- ifelse(ph[,op_code_col] %in% c(1554,1095) & !is.na(ph[,op_code_col]), TRUE, ph$cvd)
}


# Hospital records
#-----------------

if(!exists("hesin_diag")) hesin_diag <- fread("../data/d003_ukbb_hospital.tsv")
if(!exists("hesin_oper")) hesin_oper <- fread("../data/d003_ukbb_operations.tsv")


# Primary/secordary diagnosis MI
mi_diag <- hesin_diag[diag_icd9 %in% c(4109, 4129) | grepl("^I21|^I22|^I23|^I252",diag_icd10),unique(eid)]
ph$cvd <- with(ph, ifelse(iid %in% mi_diag, TRUE, ph$cvd))


# Percutaneous coronary intervention OR coronary artery bypass graft
pci_codes <- c("K49", "50", "75")
cabg_codes <- paste0("K",40:46)
ops_codes <- c(pci_codes, cabg_codes)
ops <- hesin_oper[oper3 %in% ops_codes | oper4 %in% ops_codes | posopdur %in% ops_codes | preopdur %in% ops_codes, unique(eid)]
ph$cvd <- with(ph, ifelse(iid %in% ops, TRUE, ph$cvd))



# Death registry
#---------------

# Primary or secondary cause of death is coronary artery disease

death_cause <- fread("../data/d003_ukbb_death_cause.tsv")
cad_codes <- c("I20", "I21", "I24", "I251", "I252", "I255", "I258", "I259")
cad_death <- death_cause[cause_icd10 %in% cad_codes, unique(eid)]
ph$cvd <- with(ph, ifelse(iid %in% cad_death, TRUE, ph$cvd))



#----
# Exclusions
#----

# Self-reported angina
ill_code_cols <- grep("ill_code",names(ph),value=T)
for(ill_code_col in ill_code_cols) {
  ph$cvd <- ifelse(ph$cvd==FALSE & ph[,ill_code_col]==1074 & !is.na(ph[,ill_code_col]), NA, ph$cvd)
}

# Hospital diagnosis of angina OR unstable angina
angina_diag <- hesin_diag[diag_icd9 %in% c(4139) | grepl("^I209|^I200|^I201|^I208",diag_icd10),unique(eid)]
ph$cvd <- with(ph, ifelse(ph$cvd==FALSE & iid %in% mi_diag, NA, ph$cvd))


# Self-reported angina medication
angina_med_codes <- c(
`aspirin 75mg tablet`=1140861806,
`nu-seals aspirin`=1140864860,
`aspirin`=1140868226,
`isosorbide mononitrate + aspirin`=1141164044,
`glyceryl trinitrate`=1140860834,
`gtn glyceryl trinitrate`=1140923670,
`glyceryl trinitrate patch`=1140927544,
`gtn glyceryl trinitrate patch`=1140927548,
`isosorbide mononitrate`=1140860954,
`ismn isosorbide mononitrate`=1140888762,
`ismo isosorbide mononitrate`=1140910512,
`isosorbide mononitrate product`=1141157254,
`isosorbide dinitrate`=1140861008,
`isdn isosorbide dinitrate`=1140888760,
nicorandil=1140910766)

med_code_cols <- grep("med_code",names(ph),value=T)
for(med_code_col in med_code_cols) {
  ph$cvd <- ifelse(ph$cvd==FALSE & ph[,med_code_col] %in% c(angina_med_codes) & !is.na(ph[,med_code_col]), NA, ph$cvd)
}


# Summarise

print(table(ph$cvd, useNA="ifany"))

# Write med file
fwrite(data.table(medication=names(angina_med_codes), med_code=angina_med_codes), "st003_01_ukbb_angina_med_codes.csv", sep=",", quote=FALSE, na="NA")



#----
# Other variables
#----

heading("Cholesterol medication")

medmap <- c(`Blood pressure medication`=0,
  `Cholesterol lowering medication`=1,
  `Do not know`=NA,
  `Insulin`=0,
  `None of the above`=0,
  `Prefer not to answer`=NA
)
ph$cholesterol_meds_1 <- medmap[ph$cvd_meds_1]
ph$cholesterol_meds_2 <- medmap[ph$cvd_meds_2]
ph$cholesterol_meds_3 <- medmap[ph$cvd_meds_3]

ph$cholesterol_meds <- ph$cholesterol_meds_1 | (ph$cholesterol_meds_2 & !is.na(ph$cholesterol_meds_2)) | (ph$cholesterol_meds_3 & !is.na(ph$cholesterol_meds_3))

writeLines("\nCholesterol medication:")
print(table(ph$cholesterol_meds, useNA="ifany"))
cat("\n")
print(data.table(ph)[,.N,by=c("cholesterol_meds_1","cvd_meds_1")][order(cholesterol_meds_1,-N)])
print(table(ph$cholesterol_meds_2, useNA="ifany"))
cat("\n")
print(data.table(ph)[,.N,by=c("cholesterol_meds_2","cvd_meds_2")][order(cholesterol_meds_2,-N)])
print(table(ph$cholesterol_meds_3, useNA="ifany"))
cat("\n")
print(data.table(ph)[,.N,by=c("cholesterol_meds_3","cvd_meds_3")][order(cholesterol_meds_3,-N)])



#----------------------------------------------------------------------------------------------------------------
# Survival phenotypes
#----



# Proband survival
#-----------------

heading("Proband survival")

ph_deaths <- fread("../data/d003_ukbb_deaths.tsv", select=c("eid","date_of_death"), col.names=c("iid","death_date"), data.table=FALSE)
ph_deaths <- ph_deaths[!duplicated(ph_deaths$iid),]
deceased <- ph_deaths[,"iid"]

# Add variables to main data frame
ph$dead <- ph$iid %in% deceased
ph$death_date <- ph_deaths[match(ph$iid, ph_deaths$iid), "death_date"]


# Format date variables using lubridate
ph$date_birth <- ymd(paste(floor(ph$yob),str_pad(match(ph$mob, month.name), 2, pad = "0"),"01",sep="-"))
ph$date_death <- dmy(ph$death_date)
ph$date_ass_c <- ymd(ph$date_ass_c)


# Set censoring time based on most recent death
max_date <- max(ph$date_death,na.rm=T)
end_date <- ymd("2021-11-01")
ph[which(!ph$dead), "date_death"] <- end_date


# Create study entry and exit age duration variables
ph$age_entry <- as.duration(ph$date_ass_c - ph$date_birth)
ph$age_exit <- as.duration(ph$date_death - ph$date_birth)

# Convert duration variables to years in numeric format
ph$age_entry <- as.numeric(ph$age_entry)/60/60/24/365.25
ph$age_exit <- as.numeric(ph$age_exit)/60/60/24/365.25


writeLines("\nAge entry (years):")
with(ph, print(summary(age_entry)))
writeLines("\nAge exit (years):")
with(ph, print(summary(age_exit)))
writeLines("\nCensoring time (years):")
with(ph, print(summary(age_exit-age_entry)))
cat(paste0("\nTotal person-years: ",with(ph, (format(sum(age_exit-age_entry),big.mar=","))),"\n"))




# Parent survival
#-----------------


heading("Parent survival")

survmap <- c(`No`=TRUE,
`Yes`=FALSE,
`Do not know`=NA,
`Prefer not to answer`=NA
)

# Self-reported father survival
ph$fath_dead <- survmap[ph$fath_alive]

writeLines("\nSelf-reported father deaths:")
print(table(ph$fath_dead, useNA="ifany"))
cat("\n")
print(data.table(ph)[,.N,by=c("fath_dead","fath_alive")][order(fath_dead,-N)])


# Create censored age variable
ph$fath_age <- with(ph, ifelse(is.na(fath_age), fath_age_d, fath_age))

# Exclude early deaths
writeLines(paste0("Excluding ",sum(ph$fath_age < 40, na.rm=T)," fathers aged < 40"))
ph$fath_age <- ifelse(ph$fath_age < 40, NA, ph$fath_age) 
ph$fath_age <- with(ph, tukey_qc(fath_age))

print(data.table(ph)[!is.na(fath_dead),.(.N, 
	mean_fath_age=mean(fath_age, na.rm=T),
	min_fath_age=min(fath_age, na.rm=T),
	max_fath_age=max(fath_age, na.rm=T)),by=c("fath_dead")][order(fath_dead,-N)])



# Self-reported mother survival
ph$moth_dead <- survmap[ph$moth_alive]

writeLines("\nSelf-reported mother deaths:")
print(table(ph$moth_dead, useNA="ifany"))
cat("\n")
print(data.table(ph)[,.N,by=c("moth_dead","moth_alive")][order(moth_dead,-N)])


# Create censored age variable
ph$moth_age <- with(ph, ifelse(is.na(moth_age), moth_age_d, moth_age))

# Exclude early deaths
writeLines(paste0("Excluding ",sum(ph$moth_age < 40, na.rm=T)," mothers aged < 40"))
ph$moth_age <- ifelse(ph$moth_age < 40, NA, ph$moth_age) 
ph$moth_age <- with(ph, tukey_qc(moth_age))

print(data.table(ph)[!is.na(moth_dead),.(.N, 
	mean_moth_age=mean(moth_age, na.rm=T),
	min_moth_age=min(moth_age, na.rm=T),
	max_moth_age=max(moth_age, na.rm=T)),by=c("moth_dead")][order(moth_dead,-N)])




#----------------------------------------------------------------------------------------------------------------
# Write
#----


heading("Write to file")

variables_of_interest <- c(
	"iid","bmi","age","array","cvd","age_entry","age_exit","dead",
  "sys_blood_pressure","dia_blood_pressure",
	"hypertension","n_days_moderate_exercise","avg_houshold_income",
	"smoking_history","further_education","current_employment",
	"alcohol_intake_freq",
  "fath_heart_disease","moth_heart_disease","bro_heart_disease",
  "fath_heart_disease2","moth_heart_disease2","bro_heart_disease2",
  "fath_hypertension","moth_hypertension","bro_hypertension",
  "fath_hypertension2","moth_hypertension2","bro_hypertension2",
	"fath_age","fath_dead","moth_age","moth_dead",
	"antihypertensive_meds","cholesterol_meds",
  paste0("pc",1:40))

fwrite(ph[,variables_of_interest], "st003_01_ukbb_base_data.tsv", sep="\t", quote=TRUE, na="NA")
system("gzip -9f st003_01_ukbb_base_data.tsv")


writeLines("\nExported all variables of interest to file 'st003_01_ukbb_base_data.tsv.gz'")