#!/bin/bash

REACTION=$1
EVENTS=$2
WRAPPER_VERSION=gluex_MCwrapper_srcct-v2.10.1

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml
export MCWRAPPER_CENTRAL=/work/halld2/home/boyu/src_software_builds/${WRAPPER_VERSION}
export PATH=/work/halld2/home/boyu/src_software_builds/${WRAPPER_VERSION}:${PATH}

python get_configs.py ${REACTION}             # Generate the config files for each individual run.
python get_number.py ${REACTION} ${EVENTS}    # Generate a text file of the number of events for each individual run

while read -r run number;
do
    echo ${run} ${number}
    gluex_MC.py configs/wrapper_${REACTION}_${run}.cfg ${run} ${number} per_file=1000000 batch=1 logdir=/farm_out/boyu/src_analysis/sim
done < list.txt

swif2 run -workflow src_analysis_sim_deuterium
swif2 run -workflow src_analysis_sim_helium
swif2 run -workflow src_analysis_sim_carbon

rm list.txt
rm configs/wrapper_${REACTION}_*.cfg