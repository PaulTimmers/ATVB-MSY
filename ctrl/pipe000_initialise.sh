#!/bin/bash

#----
# Setup environment
#----

MAIN=/exports/igmm/eddie/wilson-lab/projects/prj_180_ygen_ukbb
mkdir -p ${MAIN}
cd ${MAIN}


#----
# Get software
#----

mkdir -p ${MAIN}/scripts
cd ${MAIN}/scripts

if [[ ! -d ${MAIN}/scripts/yhaplo ]]; then
    git clone https://github.com/23andMe/yhaplo.git
    module load igmm/apps/python/3.7.3
    virtualenv venv
    source venv/bin/activate
    pip install --upgrade pip
    cd yhaplo
    pip install --editable .
    pip install biopython
    yhaplo -ex # Test
fi


#----
# Get data
#----

mkdir -p ${MAIN}/data/
cd ${MAIN}/data/


# ISOGG snp list
#---------------

source ../scripts/venv/bin/activate
yhaplo
mv output/isogg.snps.unique.*.txt d001_isogg_snps_hg19.tsv
rm -r output


# UKBB snp list
#--------------

ln -sf /exports/igmm/eddie/UK-BioBank-Genotype/genotypes/ukb_snp_chrY_v2.bim d002_ukbb_snps.hg19.bim
ln -sf /exports/igmm/eddie/UK-BioBank-Genotype/genotypes/ukb_cal_chrY_v2.bed d002_ukbb_snps.hg19.bed
ln -sf /exports/igmm/eddie/UK-BioBank-proj19655/genotypes/array/ukbb_proj19655.fam d002_ukbb_snps.hg19.fam


# UKBB phenotypes
#----------------

# All main variables
ln -sf /exports/igmm/eddie/UK-BioBank-proj19655/phenotypes/2021_04_27/dataset_45130/p03_usable_data/ukb45130.tsv d003_ukbb_data.tsv

# Sample qc data
ln -sf /exports/igmm/eddie/UK-BioBank-proj19655/processing/ukbb_proj19655_sample_qc.tsv d003_ukbb_sample_qc.tsv

# Relatedness data
ln -sf /exports/igmm/eddie/UK-BioBank-proj19655/processing/ukb19655_rel_s488363.dat d003_ukbb_rel.dat

# Blood biochemistry data
ln -sf /exports/igmm/eddie/UK-BioBank-proj19655/phenotypes/2019_03_27/dataset_27719/p03_usable_data/ukb27719.tsv d003_ukbb_biochem.tsv

# Death records
ln -sf /exports/igmm/eddie/UK-BioBank-proj19655/phenotypes/2021_12_10/dataportal/data/death.txt d003_ukbb_deaths.tsv
ln -sf /exports/igmm/eddie/UK-BioBank-proj19655/phenotypes/2021_12_10/dataportal/data/death_cause.txt d003_ukbb_death_cause.tsv

# Hospital records
ln -sf /exports/igmm/eddie/UK-BioBank-proj19655/phenotypes/2021_12_10/dataportal/data/hesin_diag.txt d003_ukbb_hospital.tsv

# Operation records
ln -sf /exports/igmm/eddie/UK-BioBank-proj19655/phenotypes/2021_12_10/dataportal/data/hesin_oper.txt d003_ukbb_operations.tsv


#----
# Get shapefiles
#----


mkdir -p ${MAIN}/data/d003_ukbb_geography
cd ${MAIN}/data/d003_ukbb_geography


# All ShapeFiles from InFuse UK
# http://infuse.ukdataservice.ac.uk/help/definitions/2011geographies/index.html

# UK countries
wget "https://borders.ukdataservice.ac.uk/ukborders/easy_download/prebuilt/shape/infuse_ctry_2011_clipped.zip" -O ctry_2011.zip
unzip ctry_2011.zip && rm ctry_2011.zip

# # Republic of Ireland
# wget "http://census.cso.ie/censusasp/saps/boundaries/Census2011_NUTS2_generalised20m.zip" -O ireland.zip
# unzip ireland.zip && rm ireland.zip


# English Regions
wget "https://borders.ukdataservice.ac.uk/ukborders/easy_download/prebuilt/shape/infuse_rgn_2011_clipped.zip" -O rgn_2011.zip
unzip rgn_2011.zip && rm rgn_2011.zip


# Local Authorities
wget 'https://borders.ukdataservice.ac.uk/ukborders/easy_download/prebuilt/shape/infuse_dist_lyr_2011_clipped.zip' -O dist_2011.zip
unzip dist_2011.zip && rm dist_2011.zip


# Wards and Electoral Divisions
wget 'https://borders.ukdataservice.ac.uk/ukborders/easy_download/prebuilt/shape/infuse_ward_lyr_2011_clipped.zip' -O ward_2011.zip
unzip ward_2011.zip && rm ward_2011.zip


# Lower Super Output Area 
wget "https://borders.ukdataservice.ac.uk/ukborders/easy_download/prebuilt/shape/infuse_lsoa_lyr_2011_clipped.zip" -O lsoa_2011.zip
unzip lsoa_2011.zip && rm lsoa_2011.zip


# CSV file mapping smaller areas to larger ones
wget "http://infuse.ukdataservice.ac.uk/help/definitions/2011geographies/all-areas-lookup-csv.zip" -O areas_2011.zip
unzip areas_2011.zip && rm areas_2011.zip


#----
# Get townsend data
#----

# https://www.statistics.digitalresources.jisc.ac.uk/dataset/2011-uk-townsend-deprivation-scores
## Download 'LSOA Scores' (csv, 1.97 MB) as 'townsend_lsoa_2011.csv'
