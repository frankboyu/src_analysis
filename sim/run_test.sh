#!/bin/bash

REACTION=piminus_p_2H
RUN=90213
EVENTS=10000
BATCH=0

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /work/halld2/home/boyu/src_software_builds/halld_versions_srcct/recon_srcct-2021_11-dev.xml

if [${BATCH} == 0]
then
    cd output
    python $MCWRAPPER_CENTRAL/gluex_MC.py ../configs/wrapper_${REACTION}.cfg ${RUN} ${EVENTS} per_file=1000000 batch=${BATCH}
    rm -r 90*
else
    python $MCWRAPPER_CENTRAL/gluex_MC.py configs/wrapper_${REACTION}.cfg ${RUN} ${EVENTS} per_file=1000000 batch=${BATCH} logdir=/farm_out/boyu/src_analysis/sim/
fi