#!/bin/bash

#----
# Setup environment
#----

MAIN=/exports/igmm/eddie/wilson-lab/projects/prj_180_ygen_ukbb
cd ${MAIN}

mkdir -p ${MAIN}/p003_phenotypes
cd ${MAIN}/p003_phenotypes


#----
# Create base data
#----

Rscript ../scripts/a003_01_ukbb_base_data.R |& tee st003_01_ukbb_base_data.log 



#----
# Create geographic variables
#----

Rscript ../scripts/a003_02_ukbb_geo_data.R |& tee st003_02_ukbb_geo_data.log



#----
# Create biochemistry variables
#----

Rscript ../scripts/a003_03_ukbb_biochem_data.R |& tee st003_03_ukbb_biochem_data.log



#----
# Format haplogroups
#----

Rscript ../scripts/a003_04_ukbb_haplogroups.R |& tee st003_04_ukbb_haplogroups.log