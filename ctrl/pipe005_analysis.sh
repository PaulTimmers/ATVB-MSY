#!/bin/bash

#----
# Setup environment
#----

MAIN=/exports/igmm/eddie/wilson-lab/projects/prj_180_ygen_ukbb
cd ${MAIN}

mkdir -p ${MAIN}/p005_analysis
cd ${MAIN}/p005_analysis


#----
# Eales et al. model
#----

Rscript ../scripts/a005_01_eales_model.R |& tee st005_01_eales_compare.log


#----
# Geography models
#----

# Blood pressure
Rscript ../scripts/a005_02_linear_model_bp.R

# Lipids
Rscript ../scripts/a005_03_linear_model_lipids.R

# CVD
Rscript ../scripts/a005_04_logistic_model_cvd.R

# Hypertension
Rscript ../scripts/a005_05_logistic_model_hypertension.R

# Survival
Rscript ../scripts/a005_06_coxph_model_surv.R


#----
# IJK hypertension
#----

Rscript ../scripts/a005_05_logistic_model_hypertension_ijk.R