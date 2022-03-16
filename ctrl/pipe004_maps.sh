#!/bin/bash

#----
# Setup environment
#----

MAIN=/exports/igmm/eddie/wilson-lab/projects/prj_180_ygen_ukbb
cd ${MAIN}

mkdir -p ${MAIN}/p004_maps
cd ${MAIN}/p004_maps


#----
# Compile map data
#----

Rscript ../scripts/a004_01_compile_map_data.R


#----
# Compile haplo map data
#----

Rscript ../scripts/a004_02_compile_haplo_map_data.R


#----
# Plot haplogroups
#----

# Launch array job

n_haplos=`wc -l < st004_02_map_plot_haplos.txt`
qsub -N plot_haplos -t 1-${n_haplos} -l h_vmem=8G -l h_rt=00:15:00 -j y -o log_plot_haplogroup.\$TASK_ID.log -cwd -V <<"QSUB"
truncate -s 0 log_plot_haplogroup.${SGE_TASK_ID}.log
Rscript ../scripts/a004_03_plot_haplo_map.R $SGE_TASK_ID && echo "DONE" || echo "ERROR"
QSUB


# When done, clean up log files

n_haplos=`wc -l < st004_02_map_plot_haplos.txt`
for i in `eval echo {1..$n_haplos}`
do
    grep -q "DONE" log_plot_haplogroup.${i}.log && rm log_plot_haplogroup.${i}.log
done

