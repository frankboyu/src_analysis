#!/bin/bash

REACTION=piminus_p_2H
RUN=90213
EVENTS=10000
TAG=4_genout

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /work/halld2/home/boyu/src_software_builds/halld_versions_srcct/recon_srcct-2021_11-dev.xml

mkdir -p output/${REACTION}_test/${TAG}
cd output/${REACTION}_test/${TAG}
gluex_MC.py ../../../configs/wrapper_${REACTION}.cfg ${RUN} ${EVENTS} per_file=1000000 batch=0