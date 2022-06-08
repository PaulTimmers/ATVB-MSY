#!/usr/bin/env Rscript

#----
# Setup environment
#----

set.seed(1)
options(width=200)
options(digits=6)
library(data.table)

traits <- c(
"f.eid", # iid
"f.31.0.0", # sex
"f.34.0.0", # yob
"f.52.0.0", # mob
"f.53.0.0", # date_ass_c
"f.1767.0.0", # adopted
"f.1797.0.0", # fath_alive
"f.1807.0.0", # fath_age_d
"f.2946.0.0", # fath_age
paste0("f.6177.0.",0:2), # cvd_meds
paste0("f.20002.0.",0:28), # ill_code
paste0("f.20003.0.",0:47), # med_code
paste0("f.20004.0.",0:31), # op_code
paste0("f.20107.0.",0:9), # fath_ill
"f.21000.0.0", # ethnicity
"f.22001.0.0" # genetic_sex
)

trait_names <- c(
"iid",
"sex",
"yob",
"mob",
"date_ass_c",
"adopted",
"fath_alive",
"fath_age_d",
"fath_age",
paste0("cvd_meds_",1:3),
paste0("ill_code_",1:29),
paste0("med_code_",1:48),
paste0("op_code_",1:32),
paste0("fath_ill_",1:10),
"ethnicity",
"genetic_sex"
)


#----
# Load data
#----

# Load raw
rel <- fread("../data/d003_ukbb_rel.dat")
qc <- fread("../data/d003_ukbb_sample_qc.tsv")
ph <- fread("../data/d003_ukbb_data.tsv", select=traits, col.names=trait_names)
iids <- fread("../p001_qc_variants/st001_03_iids.txt")[[1]]

# Subset data to males only
males <- qc[Inferred.Gender=="M",iid]
ph <- ph[iid %in% males]

# Get list of relationship coefficients for samples in the study
rel1 <- rel[ID1 %in% iids | ID2 %in% iids]


#----
# Infer father-son pairs
#----

# Get participant date of birth (mid-month estimate)
month.len <- c(31,28,31,30,31,30,31,31,30,31,30,31); names(month.len) <- month.name
ph[, dob := as.Date(paste(yob,match(mob, month.name),floor(month.len[mob]/2), sep="-"))]

# Get living father date of birth range
ph[is.na(fath_age_d), fath_dob_min := date_ass_c - (fath_age+0.5)*365]
ph[is.na(fath_age_d), fath_dob_max := date_ass_c - (fath_age-0.5)*365]

# Get participants with correct kinship coefficient
fath_son_candidates <- rel1[Kinship > 0.2,]

# Add date of birth of participant and their father for each pair
fath_son_candidates <- ph[,.(ID2=iid, ID2_DOB=dob, ID2_FATH_DOB_MIN=fath_dob_min, ID2_FATH_DOB_MAX=fath_dob_max)][fath_son_candidates,,on=c("ID2")]
fath_son_candidates <- ph[,.(ID1=iid, ID1_DOB=dob, ID1_FATH_DOB_MIN=fath_dob_min, ID1_FATH_DOB_MAX=fath_dob_max)][fath_son_candidates,,on=c("ID1")]
fath_son_candidates <- fath_son_candidates[!(is.na(ID1_FATH_DOB_MIN)&is.na(ID2_FATH_DOB_MIN)),]

# Subset to pairs with correct father-son date of births
fath_son_candidates[,ID1_SON:=.(difftime(ID2_DOB,ID1_FATH_DOB_MIN)>0 & difftime(ID2_DOB,ID1_FATH_DOB_MAX)<0)]
fath_son_candidates[,ID2_SON:=.(difftime(ID1_DOB,ID2_FATH_DOB_MIN)>0 & difftime(ID1_DOB,ID2_FATH_DOB_MAX)<0)]
fath_list <- rbind(fath_son_candidates[ID1_SON == TRUE,.(iid=ID1, fatid=ID2)], fath_son_candidates[ID2_SON == TRUE,.(iid=ID2, fatid=ID1)])


#----
# Get father heart disease
#----

ph <- data.frame(ph[iid %in% fath_list[,c(iid,fatid)]])

ph$fath_heart_disease <- ph$fath_heart_disease2 <- ifelse(is.na(ph$fath_ill_1), NA, 1)
fath_ill_cols <- grep("fath_ill",names(ph),value=T)

# 1: Any reported incidence is a case; 2: Exclude NAs
for(fath_ill_col in fath_ill_cols) {
ph$fath_heart_disease <- ifelse(ph[,fath_ill_col] == "Heart disease" & !is.na(ph[,fath_ill_col]), 2, ph$fath_heart_disease)
ph$fath_heart_disease <- ifelse(ph$fath_heart_disease == 1 & ph[,fath_ill_col] %in% c("Prefer not to answer (group 1)") & !is.na(ph[,fath_ill_col]), NA, ph$fath_heart_disease)

ph$fath_heart_disease2 <- ifelse(ph[,fath_ill_col] == "Heart disease" & !is.na(ph[,fath_ill_col]), 2, ph$fath_heart_disease2)
ph$fath_heart_disease2 <- ifelse(ph$fath_heart_disease2 == 1 & ph[,fath_ill_col] %in% c("Prefer not to answer (group 1)", "Do not know (group 1)") & !is.na(ph[,fath_ill_col]), NA, ph$fath_heart_disease2)
}



#----
# Define heart disease
#----

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


#----
# Perform comparison
#----

fath_list1 <- data.table(ph)[,.(fatid=iid,fath_cvd=cvd)][fath_list,,on="fatid"]
fath_list2 <- data.table(ph)[,.(iid,fath_report_cvd=fath_heart_disease)][fath_list1, on="iid"]

# Heart disease vs. CAD
fath_list_hd <- fath_list2[fath_report_cvd==2,]
prop <- fath_list_hd[fath_cvd==TRUE,.N]/fath_list_hd[,.N]

dist <- array(NA, 100000)
for(i in 1:length(dist)){
	fath_list_i <- fath_list2[sample(x=1:nrow(fath_list2), size=nrow(fath_list2), replace=TRUE),]
	fath_list_hd_i <- fath_list_i[fath_report_cvd==2,]
	prop_i <- fath_list_hd_i[fath_cvd==TRUE,.N]/fath_list_hd_i[,.N]
	dist[i] <- prop_i
}
dist_ci <- quantile(dist, prob=c(1-0.975, 0.975))

writeLines(sprintf("There are %i father-son pairs, of which %i have a father with son-reported heart disease, and %i have a father with CAD.",fath_list2[,.N],fath_list_hd[,.N],fath_list_hd[fath_cvd==TRUE,.N]))
writeLines(sprintf("Proportion of father heart diseases classified as CAD: %.1f%% (95%% CI %.1f%% to %.1f%%)",100*prop,100*dist_ci[1],100*dist_ci[2]))
