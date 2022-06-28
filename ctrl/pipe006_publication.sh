#!/bin/bash

#----
# Setup environment
#----

MAIN=/exports/igmm/eddie/wilson-lab/projects/prj_180_ygen_ukbb
cd ${MAIN}

mkdir -p ${MAIN}/p006_publication
cd ${MAIN}/p006_publication


#----
# Figure 1
#----

Rscript ../scripts/a006_01_figure1.R


#----
# Figure 2
#----

Rscript ../scripts/a006_02_figure2.R


#----
# Graphical abstract
#----

Rscript  ../scripts/a006_03_graphical_abstract.R


#----
# Supplementary Tables
#----

# The various results generated in the analyses have been copied over
# into the Microsoft Excel Document titled 'st006_03_supplementary_tables.xlsx'
# and some headers have been formatted for improved readability


#----
# Supplementary Figures
#----

# The following files were compiled in a Microsoft Word Document titled 'st006_04_supplementary_figures.docx'

# 'st005_01_eales_compare.pdf' - SF1
# 'st005_02_bp_geo_model.pdf' - SF2 
# 'st005_05_hypertension_geo_model_ijk2.pdf' - SF3


#----
# Supplementary Data
#----

cd ../p004_maps/
tar -cvzf st006_05_supplementary_data.tar.gz st004_03_*_birth.png
cd ../p006_publication
mv ../p004_maps/st006_05_supplementary_data.tar.gz ./


# The tgz archive 'st006_05_supplementary_data.zip' has been uploaded to Edinburgh DataShare
# https://doi.org/10.7488/ds/3472