#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------------------------------------
# Setup environment
#----

source("../scripts/print_script_name.R")
options(width=200)
options(digits=6)


# R packages

library(data.table)


# Variables

traits <- c(
"f.eid", # iid
"f.30690.0.0", # Total cholesterol
"f.30760.0.0", # HDL cholesterol 
"f.30780.0.0", # LDL cholesterol
"f.30870.0.0" # Triglycerides
)

trait_names <- c(
"iid",
"total_cholesterol",
"hdl_cholesterol", 
"ldl_cholesterol", 
"triglycerides" 
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

ph <- fread("../data/d003_ukbb_biochem.tsv", select=traits, col.names=trait_names, data.table=FALSE)
iids <- fread("../p001_qc_variants/st001_03_iids.txt", header=FALSE, col.names=c("iid"))

# Subset to unrelated, genomically british sample
ph <- ph[ph$iid %in% iids$iid, ]


#----
# QC
#----

# Total cholesterol
#------------------

heading("Total cholesterol")

ph$total_cholesterol <- with(ph, tukey_qc(total_cholesterol))
print(summary(ph$total_cholesterol))


# HDL cholesterol
#----------------

heading("HDL cholesterol")

ph$hdl_cholesterol <- with(ph, tukey_qc(hdl_cholesterol))
print(summary(ph$hdl_cholesterol))


# LDL cholesterol
#----------------

heading("LDL cholesterol")

ph$ldl_cholesterol <- with(ph, tukey_qc(ldl_cholesterol))
print(summary(ph$ldl_cholesterol))



# Triglyerides
#-------------

heading("Triglycerides")

ph$log_triglycerides <- with(ph, tukey_qc(log(triglycerides)))
print(summary(ph$log_triglycerides))



#----
# Write results
#----

fwrite(ph, "st003_03_ukbb_biochem_data.tsv.gz", sep="\t", quote=TRUE, na="NA", compress="gzip")
