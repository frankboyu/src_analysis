#!/bin/bash

REACTION=$1
EVENTS=$2

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /work/halld2/home/boyu/src_software_builds/halld_versions_srcct/recon_srcct-2021_11-ver01_2 _1.xml

python get_configs.py ${REACTION}             # Generate the config files for each individual run.
python get_number.py ${REACTION} ${EVENTS}    # Generate a text file of the number of events for each individual run

while read -r run number;
do
    echo ${run} ${number}
    gluex_MC.py configs/wrapper_${REACTION}_${run}.cfg ${run} ${number} per_file=1000000 batch=1  logdir=/farm_out/boyu/src_analysis/sim
done < list.txt

swif2 run -workflow src_analysis_sim

rm list.txt
rm configs/wrapper_${REACTION}_*.cfg